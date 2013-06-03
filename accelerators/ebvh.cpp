
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// accelerators/ebvh.cpp*
#include "stdafx.h"
#include "accelerators/ebvh.h"
#include "probes.h"
#include "paramset.h"
#include <stack>

// EBVHAccel Local Declarations
struct EBVHPrimitiveInfo {
    EBVHPrimitiveInfo() { }
    EBVHPrimitiveInfo(int pn, const BBox &b)
        : primitiveNumber(pn), bounds(b) {
        centroid = .5f * b.pMin + .5f * b.pMax;
    }
    int primitiveNumber;
    Point centroid;
    BBox bounds;
};


struct EBVHBuildNode {
    // EBVHBuildNode Public Methods
    EBVHBuildNode() { children[0] = children[1] = NULL; }
    void InitLeaf(uint32_t first, uint32_t n, const BBox &b) {
        firstPrimOffset = first;
        nPrimitives = n;
        bounds = b;
    }
    void InitInterior(uint32_t axis, EBVHBuildNode *c0, EBVHBuildNode *c1) {
        children[0] = c0;
        children[1] = c1;
        bounds = Union(c0->bounds, c1->bounds);
        splitAxis = axis;
        nPrimitives = 0;
    }
    BBox bounds;
    EBVHBuildNode *children[2];
    uint32_t splitAxis, firstPrimOffset, nPrimitives;
};


struct CompareToMid {
    CompareToMid(int d, float m) { dim = d; mid = m; }
    int dim;
    float mid;
    bool operator()(const EBVHPrimitiveInfo &a) const {
        return a.centroid[dim] < mid;
    }
};


struct ComparePoints {
    ComparePoints(int d) { dim = d; }
    int dim;
    bool operator()(const EBVHPrimitiveInfo &a,
                    const EBVHPrimitiveInfo &b) const {
        return a.centroid[dim] < b.centroid[dim];
    }
};


struct CompareToBucket {
    CompareToBucket(int split, int num, int d, const BBox &b)
        : centroidBounds(b)
    { splitBucket = split; nBuckets = num; dim = d; }
    bool operator()(const EBVHPrimitiveInfo &p) const;

    int splitBucket, nBuckets, dim;
    const BBox &centroidBounds;
};


bool CompareToBucket::operator()(const EBVHPrimitiveInfo &p) const {
    int b = nBuckets * ((p.centroid[dim] - centroidBounds.pMin[dim]) /
            (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
    if (b == nBuckets) b = nBuckets-1;
    Assert(b >= 0 && b < nBuckets);
    return b <= splitBucket;
}


struct LinearEBVHNode {
/***************************************************************/
/*  Change this struct to use a compressed BBox representation */
/***************************************************************/
    //BBox bounds;
    unsigned int bounds_x0:10;
    unsigned int bounds_y0:10;
    unsigned int bounds_z0:10;
    unsigned int bounds_x1:10;
    unsigned int bounds_y1:10;
    unsigned int bounds_z1:10;


    union {
        uint32_t primitivesOffset;    // leaf
        uint32_t secondChildOffset;   // interior
    };

    uint8_t nPrimitives;  // 0 -> interior node
    uint8_t axis;         // interior node: xyz
    uint8_t pad[2];       // ensure 32 byte total size
};


static inline bool IntersectP(const BBox &bounds, const Ray &ray,
        const Vector &invDir, const uint32_t dirIsNeg[3]) {
    // Check for ray intersection against $x$ and $y$ slabs
    float tmin =  (bounds[  dirIsNeg[0]].x - ray.o.x) * invDir.x;
    float tmax =  (bounds[1-dirIsNeg[0]].x - ray.o.x) * invDir.x;
    float tymin = (bounds[  dirIsNeg[1]].y - ray.o.y) * invDir.y;
    float tymax = (bounds[1-dirIsNeg[1]].y - ray.o.y) * invDir.y;
    if ((tmin > tymax) || (tymin > tmax))
        return false;
    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;

    // Check for ray intersection against $z$ slab
    float tzmin = (bounds[  dirIsNeg[2]].z - ray.o.z) * invDir.z;
    float tzmax = (bounds[1-dirIsNeg[2]].z - ray.o.z) * invDir.z;
    if ((tmin > tzmax) || (tzmin > tmax))
        return false;
    if (tzmin > tmin)
        tmin = tzmin;
    if (tzmax < tmax)
        tmax = tzmax;
    return (tmin < ray.maxt) && (tmax > ray.mint);
}



// EBVHAccel Method Definitions
EBVHAccel::EBVHAccel(const vector<Reference<Primitive> > &p,
                   uint32_t mp, const string &sm) {
    maxPrimsInNode = min(255u, mp);
    for (uint32_t i = 0; i < p.size(); ++i)
        p[i]->FullyRefine(primitives);
    if (sm == "sah")         splitMethod = SPLIT_SAH;
    else if (sm == "middle") splitMethod = SPLIT_MIDDLE;
    else if (sm == "equal")  splitMethod = SPLIT_EQUAL_COUNTS;
    else {
        Warning("EBVH split method \"%s\" unknown.  Using \"sah\".",
                sm.c_str());
        splitMethod = SPLIT_SAH;
    }

    if (primitives.size() == 0) {
        nodes = NULL;
        return;
    }
    // Build EBVH from _primitives_
    PBRT_BVH_STARTED_CONSTRUCTION(this, primitives.size());

    // Initialize _buildData_ array for primitives
    vector<EBVHPrimitiveInfo> buildData;
    buildData.reserve(primitives.size());
    for (uint32_t i = 0; i < primitives.size(); ++i) {
        BBox bbox = primitives[i]->WorldBound();
        buildData.push_back(EBVHPrimitiveInfo(i, bbox));
    }

    // Recursively build EBVH tree for primitives
    MemoryArena buildArena;
    uint32_t totalNodes = 0;
    vector<Reference<Primitive> > orderedPrims;
    orderedPrims.reserve(primitives.size());
    EBVHBuildNode *root = recursiveBuild(buildArena, buildData, 0,
                                        primitives.size(), &totalNodes,
                                        orderedPrims);
    this->sceneBoundingBox = root->bounds;
    primitives.swap(orderedPrims);
        printf("EBVH created with %d nodes for %d primitives at %d bytes/node (%.2f MB)\n", totalNodes,
             (int)primitives.size(), int(sizeof(LinearEBVHNode)), float(totalNodes * sizeof(LinearEBVHNode))/(1024.f*1024.f));

    // Compute representation of depth-first traversal of EBVH tree
    nodes = AllocAligned<LinearEBVHNode>(totalNodes);
    for (uint32_t i = 0; i < totalNodes; ++i)
        new (&nodes[i]) LinearEBVHNode;
    uint32_t offset = 0;
    flattenEBVHTree(root, &offset, this->sceneBoundingBox);
    Assert(offset == totalNodes);
    PBRT_BVH_FINISHED_CONSTRUCTION(this);
}


BBox EBVHAccel::WorldBound() const {
    //return nodes ? nodes[0].bounds : BBox();
    return nodes? this->sceneBoundingBox : BBox();
}


EBVHBuildNode *EBVHAccel::recursiveBuild(MemoryArena &buildArena,
        vector<EBVHPrimitiveInfo> &buildData, uint32_t start,
        uint32_t end, uint32_t *totalNodes,
        vector<Reference<Primitive> > &orderedPrims) {
    Assert(start != end);
    (*totalNodes)++;
    EBVHBuildNode *node = buildArena.Alloc<EBVHBuildNode>();
    // Compute bounds of all primitives in EBVH node
    BBox bbox;
    for (uint32_t i = start; i < end; ++i)
        bbox = Union(bbox, buildData[i].bounds);
    uint32_t nPrimitives = end - start;
    if (nPrimitives == 1) {
        // Create leaf _EBVHBuildNode_
        uint32_t firstPrimOffset = orderedPrims.size();
        for (uint32_t i = start; i < end; ++i) {
            uint32_t primNum = buildData[i].primitiveNumber;
            orderedPrims.push_back(primitives[primNum]);
        }
        node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
    }
    else {
        // Compute bound of primitive centroids, choose split dimension _dim_
        BBox centroidBounds;
        for (uint32_t i = start; i < end; ++i)
            centroidBounds = Union(centroidBounds, buildData[i].centroid);
        int dim = centroidBounds.MaximumExtent();

        // Partition primitives into two sets and build children
        uint32_t mid = (start + end) / 2;
        if (centroidBounds.pMax[dim] == centroidBounds.pMin[dim]) {
            // Create leaf _EBVHBuildNode_
            uint32_t firstPrimOffset = orderedPrims.size();
            for (uint32_t i = start; i < end; ++i) {
                uint32_t primNum = buildData[i].primitiveNumber;
                orderedPrims.push_back(primitives[primNum]);
            }
            node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
            return node;
        }

        // Partition primitives based on _splitMethod_
        switch (splitMethod) {
        case SPLIT_MIDDLE: {
            // Partition primitives through node's midpoint
            float pmid = .5f * (centroidBounds.pMin[dim] + centroidBounds.pMax[dim]);
            EBVHPrimitiveInfo *midPtr = std::partition(&buildData[start],
                                                      &buildData[end-1]+1,
                                                      CompareToMid(dim, pmid));
            mid = midPtr - &buildData[0];
            if (mid != start && mid != end)
                // for lots of prims with large overlapping bounding boxes, this
                // may fail to partition; in that case don't break and fall through
                // to SPLIT_EQUAL_COUNTS
                break;
        }
        case SPLIT_EQUAL_COUNTS: {
            // Partition primitives into equally-sized subsets
            mid = (start + end) / 2;
            std::nth_element(&buildData[start], &buildData[mid],
                             &buildData[end-1]+1, ComparePoints(dim));
            break;
        }
        case SPLIT_SAH: default: {
            // Partition primitives using approximate SAH
            if (nPrimitives <= 4) {
                // Partition primitives into equally-sized subsets
                mid = (start + end) / 2;
                std::nth_element(&buildData[start], &buildData[mid],
                                 &buildData[end-1]+1, ComparePoints(dim));
            }
            else {
                // Allocate _BucketInfo_ for SAH partition buckets
                const int nBuckets = 12;
                struct BucketInfo {
                    BucketInfo() { count = 0; }
                    int count;
                    BBox bounds;
                };
                BucketInfo buckets[nBuckets];

                // Initialize _BucketInfo_ for SAH partition buckets
                for (uint32_t i = start; i < end; ++i) {
                    int b = nBuckets *
                        ((buildData[i].centroid[dim] - centroidBounds.pMin[dim]) /
                         (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
                    if (b == nBuckets) b = nBuckets-1;
                    Assert(b >= 0 && b < nBuckets);
                    buckets[b].count++;
                    buckets[b].bounds = Union(buckets[b].bounds, buildData[i].bounds);
                }

                // Compute costs for splitting after each bucket
                float cost[nBuckets-1];
                for (int i = 0; i < nBuckets-1; ++i) {
                    BBox b0, b1;
                    int count0 = 0, count1 = 0;
                    for (int j = 0; j <= i; ++j) {
                        b0 = Union(b0, buckets[j].bounds);
                        count0 += buckets[j].count;
                    }
                    for (int j = i+1; j < nBuckets; ++j) {
                        b1 = Union(b1, buckets[j].bounds);
                        count1 += buckets[j].count;
                    }
                    cost[i] = .125f + (count0*b0.SurfaceArea() + count1*b1.SurfaceArea()) /
                              bbox.SurfaceArea();
                }

                // Find bucket to split at that minimizes SAH metric
                float minCost = cost[0];
                uint32_t minCostSplit = 0;
                for (int i = 1; i < nBuckets-1; ++i) {
                    if (cost[i] < minCost) {
                        minCost = cost[i];
                        minCostSplit = i;
                    }
                }

                // Either create leaf or split primitives at selected SAH bucket
                if (nPrimitives > maxPrimsInNode ||
                    minCost < nPrimitives) {
                    EBVHPrimitiveInfo *pmid = std::partition(&buildData[start],
                        &buildData[end-1]+1,
                        CompareToBucket(minCostSplit, nBuckets, dim, centroidBounds));
                    mid = pmid - &buildData[0];
                }
                
                else {
                    // Create leaf _EBVHBuildNode_
                    uint32_t firstPrimOffset = orderedPrims.size();
                    for (uint32_t i = start; i < end; ++i) {
                        uint32_t primNum = buildData[i].primitiveNumber;
                        orderedPrims.push_back(primitives[primNum]);
                    }
                    node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
                    return node;
                }
            }
            break;
        }
        }
        node->InitInterior(dim,
                           recursiveBuild(buildArena, buildData, start, mid,
                                          totalNodes, orderedPrims),
                           recursiveBuild(buildArena, buildData, mid, end,
                                          totalNodes, orderedPrims));
    }
    return node;
}

BBox EBVHAccel::compressBBox(const BBox bbox, const BBox boxToBoundOn) {
    BBox returnVal;

    float x0 = floorf( ((float)bbox.pMin.x - (float)boxToBoundOn.pMin.x) / ((float)boxToBoundOn.pMax.x - (float)boxToBoundOn.pMin.x) * ((float)pow(2, 10)-1) );
    float y0 = floorf( ((float)bbox.pMin.y - (float)boxToBoundOn.pMin.y) / ((float)boxToBoundOn.pMax.y - (float)boxToBoundOn.pMin.y) * ((float)pow(2, 10)-1) );
    float z0 = floorf( ((float)bbox.pMin.z - (float)boxToBoundOn.pMin.z) / ((float)boxToBoundOn.pMax.z - (float)boxToBoundOn.pMin.z) * ((float)pow(2, 10)-1) );

    float x1 = ceilf( ((float)bbox.pMax.x - (float)boxToBoundOn.pMin.x) / ((float)boxToBoundOn.pMax.x - (float)boxToBoundOn.pMin.x) * ((float)pow(2, 10)-1) );
    float y1 = ceilf( ((float)bbox.pMax.y - (float)boxToBoundOn.pMin.y) / ((float)boxToBoundOn.pMax.y - (float)boxToBoundOn.pMin.y) * ((float)pow(2, 10)-1) );
    float z1 = ceilf( ((float)bbox.pMax.z - (float)boxToBoundOn.pMin.z) / ((float)boxToBoundOn.pMax.z - (float)boxToBoundOn.pMin.z) * ((float)pow(2, 10)-1) );

    returnVal = BBox(Point(x0,y0,z0), Point(x1,y1,z1));

    return returnVal;
}

BBox EBVHAccel::uncompressBBox(const BBox bbox, const BBox boxToBoundOn) {
    BBox returnVal;

    float uncomp_x0 = ((float)bbox.pMin.x / ((float)pow(2,10)-1)) * ((float)boxToBoundOn.pMax.x - (float)boxToBoundOn.pMin.x) + (float)boxToBoundOn.pMin.x;
    float uncomp_y0 = ((float)bbox.pMin.y / ((float)pow(2,10)-1)) * ((float)boxToBoundOn.pMax.y - (float)boxToBoundOn.pMin.y) + (float)boxToBoundOn.pMin.y;
    float uncomp_z0 = ((float)bbox.pMin.z / ((float)pow(2,10)-1)) * ((float)boxToBoundOn.pMax.z - (float)boxToBoundOn.pMin.z) + (float)boxToBoundOn.pMin.z;

    float uncomp_x1 = ((float)bbox.pMax.x / ((float)pow(2,10)-1)) * ((float)boxToBoundOn.pMax.x - (float)boxToBoundOn.pMin.x) + (float)boxToBoundOn.pMin.x;
    float uncomp_y1 = ((float)bbox.pMax.y / ((float)pow(2,10)-1)) * ((float)boxToBoundOn.pMax.y - (float)boxToBoundOn.pMin.y) + (float)boxToBoundOn.pMin.y;
    float uncomp_z1 = ((float)bbox.pMax.z / ((float)pow(2,10)-1)) * ((float)boxToBoundOn.pMax.z - (float)boxToBoundOn.pMin.z) + (float)boxToBoundOn.pMin.z;

    returnVal = BBox(Point(uncomp_x0,uncomp_y0,uncomp_z0), Point(uncomp_x1,uncomp_y1,uncomp_z1));

    return returnVal;
}


uint32_t EBVHAccel::flattenEBVHTree(EBVHBuildNode *node, uint32_t *offset, BBox boxToBoundOn) {
    LinearEBVHNode *linearNode = &nodes[*offset];
    
    BBox compressedBox = compressBBox(node->bounds, boxToBoundOn);
    linearNode->bounds_x0 = compressedBox.pMin.x;
    linearNode->bounds_y0 = compressedBox.pMin.y;
    linearNode->bounds_z0 = compressedBox.pMin.z;
    linearNode->bounds_x1 = compressedBox.pMax.x;
    linearNode->bounds_y1 = compressedBox.pMax.y;
    linearNode->bounds_z1 = compressedBox.pMax.z;


    BBox uncompressedBox = uncompressBBox(BBox(Point(linearNode->bounds_x0, linearNode->bounds_y0, linearNode->bounds_z0), Point(linearNode->bounds_x1, linearNode->bounds_y1, linearNode->bounds_z1)), boxToBoundOn);

    uint32_t myOffset = (*offset)++;
    if (node->nPrimitives > 0) {
        Assert(!node->children[0] && !node->children[1]);
        linearNode->primitivesOffset = node->firstPrimOffset;
        linearNode->nPrimitives = node->nPrimitives;
    }
    else {
        // Creater interior flattened EBVH node
        linearNode->axis = node->splitAxis;
        linearNode->nPrimitives = 0;
        flattenEBVHTree(node->children[0], offset, BBox(Point(uncompressedBox.pMin.x,uncompressedBox.pMin.y,uncompressedBox.pMin.z), Point(uncompressedBox.pMax.x,uncompressedBox.pMax.y,uncompressedBox.pMax.z)));
        linearNode->secondChildOffset = flattenEBVHTree(node->children[1],
                                                       offset, BBox(Point(uncompressedBox.pMin.x,uncompressedBox.pMin.y,uncompressedBox.pMin.z), Point(uncompressedBox.pMax.x,uncompressedBox.pMax.y,uncompressedBox.pMax.z)));
    }
    return myOffset;
}


EBVHAccel::~EBVHAccel() {
    FreeAligned(nodes);
}


bool EBVHAccel::Intersect(const Ray &ray, Intersection *isect) const {
    if (!nodes) return false;
    PBRT_BVH_INTERSECTION_STARTED(const_cast<EBVHAccel *>(this), const_cast<Ray *>(&ray));
    bool hit = false;
    Vector invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
    uint32_t dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };
    // Follow ray through EBVH nodes to find primitive intersections
    uint32_t todoOffset = 0, nodeNum = 0;
    uint32_t todo[64];

    std::stack<BBox> stackOfBoxes;
    stackOfBoxes.push(this->sceneBoundingBox);

    while (true) {
        LinearEBVHNode *node = &nodes[nodeNum];

        BBox myBBox = uncompressBBox(BBox(Point(node->bounds_x0, node->bounds_y0, node->bounds_z0), Point(node->bounds_x1, node->bounds_y1, node->bounds_z1)), stackOfBoxes.top());

        // Check ray against EBVH node
        if (::IntersectP(myBBox, ray, invDir, dirIsNeg)) {
            if (node->nPrimitives > 0) {
                // Intersect ray with primitives in leaf EBVH node
                PBRT_BVH_INTERSECTION_TRAVERSED_LEAF_NODE(const_cast<LinearEBVHNode *>(node));
                for (uint32_t i = 0; i < node->nPrimitives; ++i)
                {
                    PBRT_BVH_INTERSECTION_PRIMITIVE_TEST(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
                    if (primitives[node->primitivesOffset+i]->Intersect(ray, isect))
                    {
                        PBRT_BVH_INTERSECTION_PRIMITIVE_HIT(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
                        hit = true;
                    }
                    else {
                        PBRT_BVH_INTERSECTION_PRIMITIVE_MISSED(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
                   }
                }
                if (todoOffset == 0) break;
                todoOffset--;
                nodeNum = todo[todoOffset];
                stackOfBoxes.pop();
            }
            else {
                // Put far EBVH node on _todo_ stack, advance to near node
                PBRT_BVH_INTERSECTION_TRAVERSED_INTERIOR_NODE(const_cast<LinearEBVHNode *>(node));
                if (dirIsNeg[node->axis]) {
                   todo[todoOffset] = nodeNum + 1;
                   nodeNum = node->secondChildOffset;
                   stackOfBoxes.pop();
                   stackOfBoxes.push(myBBox);
                   stackOfBoxes.push(myBBox);
                   todoOffset++;
                }
                else {
                   todo[todoOffset] = node->secondChildOffset;
                   nodeNum = nodeNum + 1;
                   stackOfBoxes.pop();
                   stackOfBoxes.push(myBBox);
                   stackOfBoxes.push(myBBox);
                   todoOffset++;
                }
            }
        }
        else {
            if (todoOffset == 0) break;
            todoOffset--;
            nodeNum = todo[todoOffset];
            stackOfBoxes.pop();
        }
    }
    PBRT_BVH_INTERSECTION_FINISHED();
    return hit;
}


bool EBVHAccel::IntersectP(const Ray &ray) const {
    if (!nodes) return false;
    PBRT_BVH_INTERSECTIONP_STARTED(const_cast<EBVHAccel *>(this), const_cast<Ray *>(&ray));
    Vector invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
    uint32_t dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };
    uint32_t todo[64];
    uint32_t todoOffset = 0, nodeNum = 0;

    std::stack<BBox> stackOfBoxes;
    stackOfBoxes.push(this->sceneBoundingBox);

    while (true) {
        const LinearEBVHNode *node = &nodes[nodeNum];

        BBox myBBox = uncompressBBox(BBox(Point(node->bounds_x0, node->bounds_y0, node->bounds_z0), Point(node->bounds_x1, node->bounds_y1, node->bounds_z1)), stackOfBoxes.top());

        if (::IntersectP(myBBox, ray, invDir, dirIsNeg)) {
            // Process EBVH node _node_ for traversal
            if (node->nPrimitives > 0) {
                PBRT_BVH_INTERSECTIONP_TRAVERSED_LEAF_NODE(const_cast<LinearEBVHNode *>(node));
                for (uint32_t i = 0; i < node->nPrimitives; ++i) {
                    PBRT_BVH_INTERSECTIONP_PRIMITIVE_TEST(const_cast<Primitive *>(primitives[node->primitivesOffset + i].GetPtr()));
                    if (primitives[node->primitivesOffset+i]->IntersectP(ray)) {
                        PBRT_BVH_INTERSECTIONP_PRIMITIVE_HIT(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
                        return true;
                    }
                    else {
                        PBRT_BVH_INTERSECTIONP_PRIMITIVE_MISSED(const_cast<Primitive *>(primitives[node->primitivesOffset + i].GetPtr()));
                    }
                }
                if (todoOffset == 0) break;
                todoOffset--;
                nodeNum = todo[todoOffset];
                stackOfBoxes.pop();
            }
            else {
                PBRT_BVH_INTERSECTIONP_TRAVERSED_INTERIOR_NODE(const_cast<LinearEBVHNode *>(node));
                if (dirIsNeg[node->axis]) {
                   /// second child first
                   todo[todoOffset] = nodeNum + 1;
                   nodeNum = node->secondChildOffset;

                   todoOffset++;
                   stackOfBoxes.pop();
                   stackOfBoxes.push(myBBox);
                   stackOfBoxes.push(myBBox);
                }
                else {
                   todo[todoOffset] = node->secondChildOffset;
                   nodeNum = nodeNum + 1;
                   todoOffset++;
                   stackOfBoxes.pop();
                   stackOfBoxes.push(myBBox);
                   stackOfBoxes.push(myBBox);
                }
            }
        }
        else {
            if (todoOffset == 0) break;
            todoOffset--;
            nodeNum = todo[todoOffset];
            stackOfBoxes.pop();
        }
    }
    PBRT_BVH_INTERSECTIONP_FINISHED();
    return false;
}


EBVHAccel *CreateEBVHAccelerator(const vector<Reference<Primitive> > &prims,
        const ParamSet &ps) {
    string splitMethod = ps.FindOneString("splitmethod", "sah");
    uint32_t maxPrimsInNode = ps.FindOneInt("maxnodeprims", 4);
    return new EBVHAccel(prims, maxPrimsInNode, splitMethod);
}


