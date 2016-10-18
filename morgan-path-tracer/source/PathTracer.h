#pragma once

#include <G3D/G3DAll.h>

class PathTracer : public ReferenceCountedObject {
public:

    class Options {
    public:
        
        int         raysPerPixel = 64;

        /** 1 = direct illumination */
        int         maxScatteringEvents = 5;

        Options()
#		ifdef G3D_DEBUG
			: raysPerPixel(1),
			maxScatteringEvents(1)
#		endif
		{}
    };

protected:
    typedef Point2                              PixelCoord;

    mutable TriTree					            m_triTree;
    mutable RealTime                            m_lastTreeBuildTime;

    /** For the active trace */
    mutable Options						        m_options;
    shared_ptr<Scene>                           m_scene;
    mutable shared_ptr<CubeMap>                 m_skybox;

    static const Ray                            s_degenerateRay;

    PathTracer();

    Radiance3 skyRadiance(const Vector3& direction) const;

    /** Produces a buffer of eye rays, stored in raster order in the preallocated rayBuffer. 
        \param castThroughCenter When true (for the first ray at each pixel), cast the ray through
               the pixel center to make images look less noisy.*/
	void generateEyeRays
       (int                                     width,
        int                                     height,
        const shared_ptr<Camera>&               camera,
        Array<Ray>&                             rayBuffer,
		bool									randomSubpixelPosition,
        Array<PixelCoord>&                      pixelCoordBuffer,
        const shared_ptr<Image>&                weightSumImage) const;

    /** In a properly modeled scene with area lights and no duplicating point lights, we should only count this 
        term on the first bounce. However, we're only going to sample point lights explicitly, so we need emissive
        on every bounce. Scenes like G3D cornell box where there are both point and emissives in the same location
        will get brighter than expected as a result. */
    void addEmissive
       (const Array<Ray>&                       rayFromEye,
        const Array<shared_ptr<Surfel>>&        surfelBuffer, 
        const Array<Color3>&                    modulationBuffer,
        const shared_ptr<Image>&                radianceImage,
        const Array<PixelCoord>&                pixelCoordBuffer) const;

    /** Choose what light surface to sample, storing the corresponding shadow ray and biradiance value */
    void chooseLights
       (const Array<shared_ptr<Surfel>>&        surfelBuffer,
        const Array<shared_ptr<Light>>&         lightArray,
        const Array<Ray>&                       rayBuffer,
        Array<Radiance3>&                       directBuffer,
        Array<Ray>&                             shadowRayBuffer) const;

    /** Apply the BSDF for each surfel to the biradiance in the corresponding light (unless shadowed),
        modulate as specified, and add to the image. Emissive light is only added for primary surfaces
        since it is already accounted for by explicit light sampling. */
    void shade
       (const Array<shared_ptr<Surfel>>&        surfelBuffer,
        const Array<Ray>&                       rayFromEye,
        const Array<Ray>&                       rayFromLight,
        const Array<bool>&                      lightShadowedBuffer,
        const Array<Radiance3>&                 directBuffer,
        const Array<Color3>&                    modulationBuffer,
        const shared_ptr<Image>&                radianceImage,
        const Array<PixelCoord>&                pixelCoordBuffer) const; 

    void importanceSampleLight
       (const Array<shared_ptr<Light>>&         lightArray,
        const Vector3&                          w_o,
        const shared_ptr<Surfel>&               surfel,
        Radiance3&                              weightedScatteredRadiance,
        Point3&                                 lightPosition) const;

    /** Compute the next bounce direction by mutating rayBuffer, and then multiply the modulationBuffer by each
        surfel's BSDF divided by sampling probability */
    void scatterRays
       (const Array<shared_ptr<Surfel>>&        surfelBuffer, 
        Array<Ray>&                             rayBuffer,
        Array<Color3>&                          modulationBuffer) const;

public:

    static shared_ptr<PathTracer> create();

    /** Replaces the previous scene.*/
    void setScene(const shared_ptr<Scene>& scene);

    /** Assumes that the scene has been previously set. Only rebuilds the tree
        if the scene has changed. 

        \param statusCallback Function called periodically to update the GUI with the rendering progress. Arguments are percentage (between 0 and 1) and an arbitrary message string.
      */
    void traceImage(const shared_ptr<Image>& radianceImage, const shared_ptr<Camera>& camera, const Options& options, const std::function<void(const String&, float)>& statusCallback = nullptr) const;

};
