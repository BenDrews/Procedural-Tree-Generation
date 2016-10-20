#include "PathTracer.h"
#include "bsdf.h" // For testing TODO: REMOVE

const Ray PathTracer::s_degenerateRay = Ray(Point3(1000, 1000, 1000), Vector3::unitY(), 0.0f, 0.01f);

PathTracer::PathTracer() {}

shared_ptr<PathTracer> PathTracer::create() {
    return createShared<PathTracer>();
}


void PathTracer::setScene(const shared_ptr<Scene>& scene) {
    m_scene = scene;

    // Force the scene to rebuild
    m_lastTreeBuildTime = 0;
}


Radiance3 PathTracer::skyRadiance(const Vector3& direction) const {
    if (m_skybox) {
        return m_skybox->bilinear(-direction);
    } else {
        const float r = 1.0f - powf(max(-direction.y, 0.0f), 0.8f) * 0.7f;
        return Radiance3(r, lerp(r, 1.0f, 0.6f), 1.0f) * (0.4f + max(-direction.y, 0.0f) * 0.3f);
    }
}

// Below this, paths are terminated
static const float minModulation = 0.05f;

void PathTracer::traceImage(const shared_ptr<Image>& radianceImage, const shared_ptr<Camera>& camera, const Options& options, const std::function<void(const String&, float)>& statusCallback) const {
    debugAssert(notNull(m_scene));

    if (max(m_scene->lastEditingTime(), m_scene->lastStructuralChangeTime(), m_scene->lastVisibleChangeTime()) > m_lastTreeBuildTime) {
        // Reset the tree
        m_triTree.setContents(m_scene);
        m_lastTreeBuildTime = System::time();
        m_skybox = m_scene->skyboxAsCubeMap();
    }

    // Intentionally copy the array
    Array<shared_ptr<Light>> lightArray = m_scene->lightingEnvironment().lightArray;

	// Remove disabled lights
	for (int i = lightArray.size() - 1; i >= 0; --i) {
		if (! lightArray[i]->enabled()) {
			lightArray.fastRemove(i);
		}
	}

    // Save the options to avoid passing them through every method and obscuring the algorithm
    m_options = options;

    // How much the surfaces between the eye and the current path node
    // have already modulated the contribution observed due to the BSDF.
    // Initialized based on the number of rays per pixel.
    Array< Color3 >                 modulationBuffer;

    // Primary rays
    Array< Ray >                    rayBuffer;

    // Surfels hit by primary rays (may be nullptr)
    Array< shared_ptr< Surfel > >   surfelBuffer;

    // Scattered radiance due to the selected light (which may be an emissive surface), IF
    // it is visible: (B_j * |w_j . n| * f) / p_j
    // The actual light position is implicitly encoded in the shadowRayBuffer.
    Array< Radiance3 >              directBuffer;

    // Shadow rays
    Array< Ray >                    shadowRayBuffer;

    // False if the light is visible
    Array< bool >                   lightShadowedBuffer;

    // Which pixel this data should write to. Each pixel is represented
    // at most once in the buffer. As rays are retired, this buffer
    // becomes necessary for indexing.
    Array<PixelCoord>               pixelCoordBuffer;

    // Total contribution, taking individual ray filter footprints into account
    const shared_ptr<Image>&        weightSumImage = Image::create(radianceImage->width(), radianceImage->height(), ImageFormat::R32F());

    // Resize all buffers for one sample per pixel
    const int numPixels = radianceImage->width() * radianceImage->height();
    
    // All operations act on all pixels in parallel
    radianceImage->setAll(Radiance3::zero());
    for (int r = 0; r < options.raysPerPixel; ++r) {
        modulationBuffer.resize(numPixels);
        rayBuffer.resize(numPixels);
        surfelBuffer.resize(numPixels);
        directBuffer.resize(numPixels);
        shadowRayBuffer.resize(numPixels);
        lightShadowedBuffer.resize(numPixels);
        pixelCoordBuffer.resize(numPixels);
    
        generateEyeRays(radianceImage->width(), radianceImage->height(), camera, rayBuffer, options.raysPerPixel > 1, pixelCoordBuffer, weightSumImage);
        modulationBuffer.setAll(Color3::one());
        
        // Visualize eye rays
        // for (Point2int32 P(0, 0); P.y < radianceImage->height(); ++P.y) for (P.x = 0; P.x < radianceImage->width(); ++P.x) radianceImage->set(P, Radiance3(rayBuffer[P.x + P.y * radianceImage->width()].direction() * 0.5f + Vector3::one() * 0.5f)); return;

        for (int scatteringEvents = 0; (scatteringEvents < options.maxScatteringEvents) && (surfelBuffer.size() > 0); ++scatteringEvents) {
            m_triTree.intersectRays(rayBuffer, surfelBuffer, (scatteringEvents == 0) ? TriTree::COHERENT_RAY_HINT : 0);

            addEmissive(rayBuffer, surfelBuffer, modulationBuffer, radianceImage, pixelCoordBuffer);

            // Compact buffers by removing paths that terminated (missed the entire scene)
            // This must be done serially.
            for (int i = 0; i < surfelBuffer.size(); ++i) {
                if (isNull(surfelBuffer[i]) || (modulationBuffer[i].sum() < minModulation)) {
                    modulationBuffer.fastRemove(i);
                    rayBuffer.fastRemove(i);
                    surfelBuffer.fastRemove(i);
                    directBuffer.fastRemove(i);
                    shadowRayBuffer.fastRemove(i);
                    lightShadowedBuffer.fastRemove(i);
                    pixelCoordBuffer.fastRemove(i);
                    --i;
                }
            } // for i

            // Direct lighting
            if (lightArray.size() > 0) {
                chooseLights(surfelBuffer, lightArray, rayBuffer, directBuffer, shadowRayBuffer);
                m_triTree.intersectRays(shadowRayBuffer, lightShadowedBuffer, TriTree::COHERENT_RAY_HINT | TriTree::DO_NOT_CULL_BACKFACES | TriTree::OCCLUSION_TEST_ONLY);
                //lightShadowedBuffer.setAll(false);
                shade(surfelBuffer, rayBuffer, shadowRayBuffer, lightShadowedBuffer, directBuffer, modulationBuffer, radianceImage, pixelCoordBuffer);
            }

            // Indirect lighting
            if (scatteringEvents < options.maxScatteringEvents - 1) {
                scatterRays(surfelBuffer, rayBuffer, modulationBuffer);
            }
        } // for scattering events

        if (statusCallback) { statusCallback(format("%d/%d rays/pixel", r, options.raysPerPixel), float(r) / float(options.raysPerPixel)); }
    } // for rays per pixel

    // Normalize by the weight per pixel
    Thread::runConcurrently(Point2int32(0, 0), Point2int32(radianceImage->width(), radianceImage->height()), [&](Point2int32 pix) {
        radianceImage->set(pix, radianceImage->get<Radiance3>(pix) / max(0.00001f, weightSumImage->get<Color1>(pix).value));
    });
}


void PathTracer::generateEyeRays
(int                                 width, 
 int                                 height,
 const shared_ptr<Camera>&           camera,
 Array<Ray>&                         rayBuffer,
 bool                                randomSubpixelPosition,
 Array<PixelCoord>&                  pixelCoordBuffer,
 const shared_ptr<Image>&            weightSumImage) const {

    const Rect2D viewport = Rect2D::xywh(0.0f, 0.0f, float(width), float(height));

    Thread::runConcurrently(Point2int32(0, 0), Point2int32(width, height), [&](Point2int32 point) {
        Vector2 offset(0.5f, 0.5f);
        if (randomSubpixelPosition) {
            Random& rng = Random::threadCommon();
            offset.x = rng.uniform(); offset.y = rng.uniform();
        }
        const int i = point.x + point.y * width;
	    rayBuffer[i] = camera->worldRay(float(point.x) + offset.x, float(point.y) + offset.y, viewport);

        // Camera coords put integers at top left, but image coords put them at pixel centers
        const PixelCoord& pixelCoord = Point2(point) + offset - Point2(0.5f, 0.5f);
        pixelCoordBuffer[i] = pixelCoord;
        weightSumImage->bilinearIncrement(pixelCoord, Color1(1.0f));
    });
}


void PathTracer::addEmissive
(const Array<Ray>&                   rayFromEye,
 const Array<shared_ptr<Surfel>>&    surfelBuffer, 
 const Array<Color3>&                modulationBuffer,
 const shared_ptr<Image>&            radianceImage,
 const Array<PixelCoord>&            pixelCoordBuffer) const {
    
    const int width  = radianceImage->width(), height = radianceImage->height();
    Thread::runConcurrently(0, rayFromEye.length(), [&](int i) {
        const shared_ptr<Surfel>& surfel = surfelBuffer[i];
            
        const Vector3& w_o = -rayFromEye[i].direction();
        const Radiance3& L_e = notNull(surfel) ? surfel->emittedRadiance(w_o) : skyRadiance(w_o);

        if (L_e.nonZero()) {
            // Don't incur the memory transaction cost for the common case of no emissive
            radianceImage->bilinearIncrement(pixelCoordBuffer[i], L_e * modulationBuffer[i]);
        }
    });
}


void PathTracer::importanceSampleLight
(const Array<shared_ptr<Light>>&             lightArray,
 const Vector3&                              w_o,
 const shared_ptr<Surfel>&                   surfel,
 Radiance3&                                  weightedScatteredRadiance,
 Point3&                                     lightPosition) const {

    const Point3&  X = surfel->position;
    const Vector3& n = surfel->shadingNormal;

    if (lightArray.size() == 1) {
        // There is only one light, so of course we will sample it
        const shared_ptr<Light>& light = lightArray[0];
        lightPosition = light->position().xyz();
        const Vector3& w_i = (lightPosition - X).direction();    
        weightedScatteredRadiance = light->biradiance(X) * surfel->finiteScatteringDensity(w_i, w_o) * abs(w_i.dot(n));
    } else {

        // Compute the biradiance for each light. Assume a small number of lights (e.g., 10)
        SmallArray<Radiance3, 10> radianceForLight;
        radianceForLight.resize(lightArray.size());
        Radiance totalRadiance = 0.0f;

        for (int j = 0; j < lightArray.size(); ++j) {
            const shared_ptr<Light>& light = lightArray[j];
            const Vector3& w_i = (light->position().xyz() - X).direction();    
            const Radiance3& L = light->biradiance(X) * surfel->finiteScatteringDensity(w_i, w_o) * abs(w_i.dot(n));
            radianceForLight[j] = L;
            totalRadiance += L.sum();
        }

        // Choose light index j proportional to the light's relative biradiance. Note that 
        // we always select the last light if we slightly overshot due to roundoff. In scenes
        // with only one light, we always choose that light of course.
        int j = 0;
        Random& rng = Random::threadCommon();
        for (float r = rng.uniform(0, totalRadiance); j < lightArray.size() - 1; ++j) {
            r -= radianceForLight[j].sum();
            if (r <= 0.0f) { break; }
        }

        // Choose proportional to total biradiance
        const Radiance3& L = radianceForLight[min(j, radianceForLight.size())];
        // Net radiance = L / (probability of choosing this light)
        weightedScatteredRadiance = L * (totalRadiance / max(L.sum(), 0.0001f));
        lightPosition = lightArray[j]->position().xyz();
    }
}


void PathTracer::chooseLights
(const Array<shared_ptr<Surfel>>&    surfelBuffer, 
 const Array<shared_ptr<Light>>&     lightArray,
 const Array<Ray>&                   rayBuffer,
 Array<Radiance3>&                   directBuffer,
 Array<Ray>&                         shadowRayBuffer) const {

    const float epsilon = 1e-3f;

    Thread::runConcurrently(0, surfelBuffer.size(), [&](int i) {
        const shared_ptr<Surfel>& surfel = surfelBuffer[i];

        debugAssert(notNull(surfel));
        Point3 lightPosition;
        Radiance3& L_sd = directBuffer[i];
        importanceSampleLight(lightArray, -rayBuffer[i].direction(), surfel, L_sd, lightPosition);

        // Cast shadow rays from the light to the surface for more coherence in scenes
        // with few lights (i.e., where many pixels are casting from the same lights)
        if (L_sd.nonZero()) {
            const Vector3& delta = surfel->position + surfel->geometricNormal * epsilon - lightPosition;
            const float distance = delta.length();
            shadowRayBuffer[i] = Ray::fromOriginAndDirection(lightPosition, delta * (1.0f / distance), 0.0f, distance - epsilon);
        } else {
            // Backface: create a degenerate ray
            shadowRayBuffer[i] = s_degenerateRay;
            directBuffer[i] = Radiance3::zero();
        }
    });
}


void PathTracer::shade
(const Array<shared_ptr<Surfel>>&        surfelBuffer,
 const Array<Ray>&                       rayFromEye,
 const Array<Ray>&                       rayFromLight,
 const Array<bool>&                      lightShadowedBuffer,
 const Array<Radiance3>&                 directBuffer,
 const Array<Color3>&                    modulationBuffer,
 const shared_ptr<Image>&                radianceImage,
 const Array<PixelCoord>&                pixelCoordBuffer) const {

    Thread::runConcurrently(0, surfelBuffer.size(), [&](int i) {
        const shared_ptr<Surfel>& surfel = surfelBuffer[i];
        
        debugAssertM(notNull(surfel), "Null surfels should have been compacted before shading"); 

        const Radiance3& L_sd = directBuffer[i];
        if (L_sd.nonZero() && ! lightShadowedBuffer[i]) {
            // Gamma correct textures!
            Color3& c = dynamic_pointer_cast<UniversalSurfel>(surfel)->lambertianReflectivity;
            c *= c;

#			if 1 // Regular implementation
			    const Radiance3& L = L_sd;
			    debugAssert(L.min() >= 0.0f);
#			elif 0 // Visualize lambertian term
                const Vector3& w_i = -rayFromLight[i].direction();
                const Vector3& w_o = -rayFromEye[i].direction();
                const Vector3& n   = surfel->shadingNormal;
                const shared_ptr<UniversalSurfel>& us = dynamic_pointer_cast<UniversalSurfel>(surfel);
			    const Radiance3& L = us->lambertianReflectivity / pif();
#			elif 0 // Visualize hit points
			    const Radiance3& L = Radiance3(surfel->position * 0.3f + Point3::one() * 0.5f));
#			else // Visualize normals
			    const Radiance3& L = Radiance3(surfel->geometricNormal * 0.5f + Point3::one() * 0.5f));
#			endif

			radianceImage->bilinearIncrement(pixelCoordBuffer[i], L * modulationBuffer[i]);	
        }
    });
}


void PathTracer::scatterRays
(const Array<shared_ptr<Surfel>>&        surfelBuffer, 
 Array<Ray>&                             rayBuffer,
 Array<Color3>&                          modulationBuffer) const {

    static const float epsilon = 1e-4f;

    Thread::runConcurrently(0, surfelBuffer.size(), [&](int i) {
        const shared_ptr<Surfel>& surfel = surfelBuffer[i];
        debugAssertM(notNull(surfel), "Null surfels should have been compacted before scattering"); 

        // Direction that the light went OUT, eventually towards the eye
        const Vector3& w_o = -rayBuffer[i].direction();

		// Let p(w) be the PDF we're actually sampling to produce an output vector
		// Let g(w) = f(w, w') |w . n|   [Note: g() is *not* a PDF; just the arbitrary integrand]
		//
        // w = sample with respect to p(w)
 		// weight = g(w) / p(w)
        Color3 weight;

        // Direction that light came in, being sampled
        Vector3 w_i;

#		if 1 // G3D scattering
			surfel->scatter(PathDirection::EYE_TO_SOURCE, w_o, false, Random::threadCommon(), weight, w_i);
#		else // Explicit BSDF for debugging BSDFs. This has nothing to do with the path tracer assignment!
            // scatterDBRDF
            //scatterDisney
            // scatterPeteCone
            // scatterBlinnPhong
            // scatterHackedBlinnPhong
            SimpleBSDF::scatter
            (dynamic_pointer_cast<UniversalSurfel>(surfel), w_o, Random::threadCommon(), w_i, weight);
#        endif


        modulationBuffer[i] *= weight;

        if ((modulationBuffer[i].sum() < minModulation) || w_i.isNaN()) {
            // This ray didn't scatter; it will be culled after the next pass, so avoid
            // the cost of a real ray cast on it
            rayBuffer[i] = s_degenerateRay;
        } else {
            rayBuffer[i] = Ray::fromOriginAndDirection(surfel->position + surfel->geometricNormal * epsilon * sign(w_i.dot(surfel->geometricNormal)), w_i);
        }
    });
}
 
