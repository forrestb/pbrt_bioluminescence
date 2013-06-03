
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

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CORE_API_H
#define PBRT_CORE_API_H

// core/api.h*
#include "pbrt.h"

#include "stdafx.h"
#include "api.h"
#include "parallel.h"
#include "paramset.h"
#include "spectrum.h"
#include "scene.h"
#include "renderer.h"
#include "film.h"
#include "volume.h"
#include "probes.h"

// API Additional Headers
#include "accelerators/bvh.h"
#include "accelerators/ebvh.h" 
#include "accelerators/grid.h"
#include "accelerators/kdtreeaccel.h"
#include "cameras/environment.h"
#include "cameras/orthographic.h"
#include "cameras/perspective.h"
#include "film/image.h"
#include "filters/box.h"
#include "filters/gaussian.h"
#include "filters/mitchell.h"
#include "filters/sinc.h"
#include "filters/triangle.h"
#include "integrators/ambientocclusion.h"
#include "integrators/diffuseprt.h"
#include "integrators/dipolesubsurface.h"
#include "integrators/directlighting.h"
#include "integrators/emission.h"
#include "integrators/glossyprt.h"
#include "integrators/igi.h"
#include "integrators/irradiancecache.h"
#include "integrators/path.h"
#include "integrators/photonmap.h"
#include "integrators/single.h"
#include "integrators/useprobes.h"
#include "integrators/whitted.h"
#include "lights/diffuse.h"
#include "lights/distant.h"
#include "lights/goniometric.h"
#include "lights/infinite.h"
#include "lights/point.h"
#include "lights/projection.h"
#include "lights/spot.h"
#include "materials/glass.h"
#include "materials/kdsubsurface.h"
#include "materials/matte.h"
#include "materials/measured.h"
#include "materials/metal.h"
#include "materials/mirror.h"
#include "materials/mixmat.h"
#include "materials/plastic.h"
#include "materials/substrate.h"
#include "materials/subsurface.h"
#include "materials/shinymetal.h"
#include "materials/translucent.h"
#include "materials/uber.h"
#include "renderers/aggregatetest.h"
#include "renderers/createprobes.h"
#include "renderers/metropolis.h"
#include "renderers/samplerrenderer.h"
#include "renderers/surfacepoints.h"
#include "samplers/adaptive.h"
#include "samplers/bestcandidate.h"
#include "samplers/halton.h"
#include "samplers/lowdiscrepancy.h"
#include "samplers/random.h"
#include "samplers/stratified.h"
#include "shapes/cone.h"
#include "shapes/cylinder.h"
#include "shapes/disk.h"
#include "shapes/heightfield.h"
#include "shapes/hyperboloid.h"
#include "shapes/loopsubdiv.h"
#include "shapes/nurbs.h"
#include "shapes/paraboloid.h"
#include "shapes/sphere.h"
#include "shapes/trianglemesh.h"
#include "textures/bilerp.h"
#include "textures/checkerboard.h"
#include "textures/constant.h"
#include "textures/dots.h"
#include "textures/fbm.h"
#include "textures/imagemap.h"
#include "textures/marble.h"
#include "textures/mix.h"
#include "textures/scale.h"
#include "textures/uv.h"
#include "textures/windy.h"
#include "textures/wrinkled.h"
#include "volumes/exponential.h"
#include "volumes/homogeneous.h"
#include "volumes/volumegrid.h"





// API Local Classes
#define MAX_TRANSFORMS 2
#define START_TRANSFORM_BITS (1 << 0)
#define END_TRANSFORM_BITS   (1 << 1)
#define ALL_TRANSFORMS_BITS  ((1 << MAX_TRANSFORMS) - 1)
struct TransformSet {
   // TransformSet Public Methods
   Transform &operator[](int i) {
       Assert(i >= 0 && i < MAX_TRANSFORMS);
       return t[i];
   }
   const Transform &operator[](int i) const { Assert(i >= 0 && i < MAX_TRANSFORMS); return t[i]; }
   friend TransformSet Inverse(const TransformSet &ts) {
       TransformSet t2;
       for (int i = 0; i < MAX_TRANSFORMS; ++i)
           t2.t[i] = Inverse(ts.t[i]);
       return t2;
   }
   bool IsAnimated() const {
       for (int i = 0; i < MAX_TRANSFORMS-1; ++i)
           if (t[i] != t[i+1]) return true;
       return false;
   }
private:
    Transform t[MAX_TRANSFORMS];
};


struct RenderOptions {
    // RenderOptions Public Methods
    RenderOptions();
    Scene *MakeScene();
    Camera *MakeCamera() const;
    Renderer *MakeRenderer() const;

    // RenderOptions Public Data
    float transformStartTime, transformEndTime;
    string FilterName;
    ParamSet FilterParams;
    string FilmName;
    ParamSet FilmParams;
    string SamplerName;
    ParamSet SamplerParams;
    string AcceleratorName;
    ParamSet AcceleratorParams;
    string RendererName;
    string SurfIntegratorName, VolIntegratorName;
    ParamSet RendererParams;
    ParamSet SurfIntegratorParams, VolIntegratorParams;
    string CameraName;
    ParamSet CameraParams;
    TransformSet CameraToWorld;
    vector<Light *> lights;
    vector<Reference<Primitive> > primitives;
    mutable vector<VolumeRegion *> volumeRegions;
    map<string, vector<Reference<Primitive> > > instances;
    vector<Reference<Primitive> > *currentInstance;
};


// API Function Declarations
void pbrtInit(const Options &opt);
void pbrtCleanup();
void pbrtIdentity();
void pbrtTranslate(float dx, float dy, float dz);
void pbrtRotate(float angle, float ax, float ay, float az);
void pbrtScale(float sx, float sy, float sz);
void pbrtLookAt(float ex, float ey, float ez,
                float lx, float ly, float lz,
                float ux, float uy, float uz);
void pbrtConcatTransform(float transform[16]);
void pbrtTransform(float transform[16]);
void pbrtCoordinateSystem(const string &);
void pbrtCoordSysTransform(const string &);
void pbrtActiveTransformAll();
void pbrtActiveTransformEndTime();
void pbrtActiveTransformStartTime();
void pbrtTransformTimes(float start, float end);
void pbrtPixelFilter(const string &name, const ParamSet &params);
void pbrtFilm(const string &type, const ParamSet &params);
void pbrtSampler(const string &name, const ParamSet &params);
void pbrtAccelerator(const string &name, const ParamSet &params);
void pbrtSurfaceIntegrator(const string &name, const ParamSet &params);
void pbrtVolumeIntegrator(const string &name, const ParamSet &params);
void pbrtRenderer(const string &name, const ParamSet &params);
void pbrtCamera(const string &, const ParamSet &cameraParams);
void pbrtWorldBegin();
void pbrtAttributeBegin();
void pbrtAttributeEnd();
void pbrtTransformBegin();
void pbrtTransformEnd();
void pbrtTexture(const string &name, const string &type,
    const string &texname, const ParamSet &params);
void pbrtMaterial(const string &name, const ParamSet &params);
void pbrtMakeNamedMaterial(const string &name, const ParamSet &params);
void pbrtNamedMaterial(const string &name);
void pbrtLightSource(const string &name, const ParamSet &params);
void pbrtAreaLightSource(const string &name, const ParamSet &params);
void pbrtShape(const string &name, const ParamSet &params);
void pbrtReverseOrientation();
void pbrtVolume(const string &name, const ParamSet &params);
void pbrtObjectBegin(const string &name);
void pbrtObjectEnd();
void pbrtObjectInstance(const string &name);
void pbrtWorldEnd();

static RenderOptions *renderOptions = NULL;

struct RenderOptions *getRenderOptions();

#endif // PBRT_CORE_API_H
