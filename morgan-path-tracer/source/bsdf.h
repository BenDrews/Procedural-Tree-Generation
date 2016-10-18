#pragma once

//////////////////////////////////////////////////////////////////////////////////////
// Simple BRDF designed to test easy implementation of impulses + arbitrary lobes

/** Reciprocal and non-energy increasing. Supports specular (a.k.a. mirror, impulse) reflection and
    refraction, glossy reflection, and lambertian reflection simultaneously. Efficient importance
    sampling and direct evaluation. Uses existing G3D material parameters. */
class SimpleBSDF {
public:
    static Color3 finiteScatteringDensity(const shared_ptr<UniversalSurfel>& surfel, const Vector3& w_o, const Vector3& w_i) {
        // Fresnel reflection at normal incidence
        const Color3& F_0 = surfel->glossyReflectionCoefficient;

        // Lambertian reflectivity (conditioned on not glossy reflected)
        const Color3& p_L = surfel->lambertianReflectivity;

        // Surface normal
        const Vector3& n = surfel->shadingNormal;

        // Half vector
        const Vector3& w_h = (w_i + w_o).directionOrZero();

        const float smoothness = surfel->smoothness;

        // Frensel reflection coefficient for this angle. Ignore fresnel
        // on surfaces that are magically set to zero reflectance.
        const Color3& F =
            F_0.nonZero() ?
            UniversalBSDF::schlickFresnel(F_0, max(0.0f, w_h.dot(w_i)), smoothness) : Color3::zero();

        // Lambertian term
        Color3 result = (Color3::one() - F) * p_L / pif();

        // Ignore mirror impulse's contribution, which is handled in getImpulses().
        if (smoothness != 1.0f) {
            // Normalized Blinn-Phong lobe
            const float m = UniversalBSDF::smoothnessToBlinnPhongExponent(smoothness);
            const float glossyLobe = pow(max(w_h.dot(n), 0.0f), m) *
                (8.0f + m) / (8.0f * pif() * square(max(w_i.dot(n), w_o.dot(n))));
            result += F * glossyLobe;
        }

        return result;
    }


    static void getImpulses(const shared_ptr<UniversalSurfel>& surfel, const Vector3& w_o, Surfel::ImpulseArray& impulseArray) {
        // Fresnel reflection at normal incidence
        const Color3& F_0 = surfel->glossyReflectionCoefficient;

        // Lambertian reflectivity (conditioned on not glossy reflected)
        const Color3& p_L = surfel->lambertianReflectivity;

        // Transmission (conditioned on not glossy or lambertian reflected)
        const Color3& T = surfel->transmissionCoefficient;

        // Surface normal
        const Vector3& n = surfel->shadingNormal;

        // The half-vector IS the normal for mirror reflection purposes
        const float smoothness = surfel->smoothness;

        // Frensel reflection coefficient for this angle. Ignore fresnel
        // on surfaces that are magically set to zero reflectance.
        const Color3& F = F_0.nonZero() ? UniversalBSDF::schlickFresnel(F_0, max(0.0f, n.dot(w_o)), smoothness) : Color3::zero();

        // Mirror reflection
        if ((smoothness == 1.0f) && F_0.nonZero()) {
            Surfel::Impulse& impulse = impulseArray.next();
            impulse.direction = w_o.reflectAbout(n);
            impulse.magnitude = F;
        }

        // Transmission
        const Color3& transmissionMagnitude = T * (Color3::one() - F) * (Color3::one() - (Color3::one() - F) * p_L);
        if (transmissionMagnitude.nonZero()) {
            const Vector3& transmissionDirection = (-w_o).refractionDirection(n, surfel->etaNeg, surfel->etaPos);

            // Test for total internal reflection
            if (transmissionDirection.nonZero()) {
                Surfel::Impulse& impulse = impulseArray.next();
                impulse.direction = transmissionDirection;
                impulse.magnitude = transmissionMagnitude;
            }
        }
    }


    static void scatter(const shared_ptr<UniversalSurfel>& surfel, const Vector3& w_o, Random& rng, Vector3& w_i, Color3& weight) {
        Surfel::ImpulseArray impulseArray;
        getImpulses(surfel, w_o, impulseArray);

        float impulseMagnitudeSum = 0.0f;
        float r = rng.uniform();
        for (int i = 0; i < impulseArray.size(); ++i) {
            const Surfel::Impulse& impulse = impulseArray[i];
            const float probabilityOfThisImpulse = impulse.magnitude.average();
            r -= probabilityOfThisImpulse;
            impulseMagnitudeSum += probabilityOfThisImpulse;
            if (r <= 0.0f) {
                w_i    = impulse.direction;
                weight = impulse.magnitude / probabilityOfThisImpulse;
                return;
            }
        }

        // Surface normal
        const Vector3& n = surfel->shadingNormal;

        // Fresnel reflection at normal incidence
        const Color3& F_0 = surfel->glossyReflectionCoefficient;

        // Estimate the fresnel term coarsely, assuming mirror reflection. This is only used
        // for estimating the relativeGlossyProbability for the pdf; error will only lead to
        // noise, not bias in the result.
        const Color3& F = F_0.nonZero() ? UniversalBSDF::schlickFresnel(F_0, max(0.0f, n.dot(w_o)), surfel->smoothness) : Color3::zero();

        // Lambertian reflectivity (conditioned on not glossy reflected)
        const Color3& p_L = surfel->lambertianReflectivity;

        // Exponent for the cosine power lobe in the PDF that we're sampling. Rolling off
        // slightly from pure Blinn-Phong appears to give faster convergence.
        const float m = UniversalBSDF::smoothnessToBlinnPhongExponent(surfel->smoothness * 0.8f);

        float relativeGlossyProbability = F_0.nonZero() ? F.average() / (F + (Color3::one() - F) * p_L).average() : 0.0f;
        float pdfValue;
        Vector3::cosHemiPlusCosPowHemiHemiRandom(w_o.reflectAbout(n), surfel->shadingNormal, m, relativeGlossyProbability, rng, w_i, pdfValue);

        // Took this branch with probability (1 - impulseMagnitudeSum)
        pdfValue *= 1.0f - impulseMagnitudeSum;
        if (pdfValue > 0.0f) {
            weight = finiteScatteringDensity(surfel, w_o, w_i) * max(w_i.dot(n), 0.0f) / pdfValue;
        } else {
            weight = Color3::zero();
        }
    }
};


//////////////////////////////////////////////////////////////////////////////////////
// Disney BRDF

float SchlickFresnel(float u) {
    // pow(m,5)
    float m = clamp(1.0f - u, 0.0f, 1.0f);
    return square(square(m)) * m;
}

float GTR1(float NdotH, float a) {
    if (a >= 1.0f) { return 1.0f / pif(); }
    float a2 = square(a);
    float t = 1.0f + (a2 - 1.0f) * NdotH * NdotH;
    return (a2 - 1.0f) / (pif() * log(a2) * t);
}

float GTR2(float NdotH, float a) {
    float a2 = square(a);
    float t = 1.0f + (a2 - 1.0f) * NdotH * NdotH;
    return a2 / (pif() * square(t));
}

float GTR2_aniso(float NdotH, float HdotX, float HdotY, float ax, float ay) {
    return 1.0f / ( pif() * ax * ay * square(square(HdotX / ax) + square(HdotY / ay) + NdotH * NdotH ));
}

float SmithG_GGX(float Ndotv, float alphaG) {
    float a = alphaG * alphaG;
    float b = Ndotv * Ndotv;
    return 1.0f / (Ndotv + sqrt(a + b - a * b));
}

// See http://blog.selfshadow.com/publications/s2012-shading-course/burley/s2012_pbs_disney_brdf_notes_v3.pdf
// and https://github.com/wdas/brdf/blob/master/src/brdfs/disney.brdf
// for documentation of material parameters. Unlike their reference code, our baseColor is in linear
// space (not gamma encoded)
//
// L is the unit vector to the light source (omega_in) in world space
// N is the unit normal in world space
// V is the vector to the eye (omega_out) in world space
// X and Y are the tangent directions in world space
//
// Nonstandard variable names are because this is a straight port of the code from Burley.
Color3 evaluateDisneyBRDF
   (Color3  baseColor,
    float   metallic,
    float   subsurface,
    float   specular,
    float   roughness,
    float   specularTint, 
    float   anisotropic,
    float   sheen,
    float   sheenTint,
    float   clearcoat, 
    float   clearcoatGloss,
    Vector3    L,
    Vector3    V,
    Vector3    N,
    Vector3    X,
    Vector3    Y) {

    float NdotL = dot(N, L);
    float NdotV = dot(N, V);
    if (NdotL < 0 || NdotV < 0) return Color3(0.0f);

    Vector3 H = normalize(L + V);
    float NdotH = dot(N, H);
    float LdotH = dot(L, H);

    float luminance = baseColor.dot(Color3(0.3f, 0.6f, 0.1f));

    // normalize luminance to isolate hue and saturation components
    Color3 Ctint = (luminance > 0.0f) ? baseColor / luminance : Color3(1.0); 
    Color3 Cspec0 = (specular * 0.08f * Color3(1.0).lerp(Ctint, specularTint)).lerp(baseColor, metallic);
    Color3 Csheen = Color3(1.0f).lerp(Ctint, sheenTint);

    // Diffuse fresnel - go from 1 at normal incidence to .5 at grazing
    // and mix in diffuse retro-reflection based on roughness
    float FL = SchlickFresnel(NdotL), FV = SchlickFresnel(NdotV);
    float Fd90 = 0.5f + 2.0f * LdotH * LdotH * roughness;
    float Fd = lerp(1.0f, Fd90, FL) * lerp(1.0f, Fd90, FV);

    // Based on Hanrahan-Krueger BRDF approximation of isotropic BSSRDF
    // 1.25 scale is used to (roughly) preserve albedo
    // Fss90 used to "flatten" retroreflection based on roughness
    float Fss90 = LdotH * LdotH * roughness;
    float Fss = lerp(1.0f, Fss90, FL) * lerp(1.0f, Fss90, FV);
    float ss = 1.25f * (Fss * (1 / (NdotL + NdotV) - 0.5f) + 0.5f);

    // Specular
    float aspect = sqrt(1.0f - anisotropic * 0.9f);
    float ax = max(0.001f, square(roughness) / aspect);
    float ay = max(0.001f, square(roughness) * aspect);
    float Ds = GTR2_aniso(NdotH, dot(H, X), dot(H, Y), ax, ay);
    float FH = SchlickFresnel(LdotH);
    Color3 Fs = Cspec0.lerp(Color3(1.0f), FH);
    float roughg = square(roughness * 0.5f + 0.5f);
    float Gs = SmithG_GGX(NdotL, roughg) * SmithG_GGX(NdotV, roughg);

    // sheen
    Color3 Fsheen = FH * sheen * Csheen;

    // clearcoat (ior = 1.5 -> F0 = 0.04)
    float Dr = GTR1(NdotH, lerp(0.1f, 0.001f, clearcoatGloss));
    float Fr = lerp(0.04f, 1.0f, FH);
    float Gr = SmithG_GGX(NdotL, 0.25f) * SmithG_GGX(NdotV, 0.25f);

    return ((1.0f / pif()) * lerp(Fd, ss, subsurface) * baseColor + Fsheen) * (1.0f - metallic) + 
        Gs * Fs * Ds + Color3(0.25f * clearcoat * Gr * Fr * Dr);
}


void scatterDisney(const shared_ptr<UniversalSurfel>& surfel, const Vector3& w_o, Random& rng, Vector3& w_i, Color3& weight) {
    const Color3& p_L = surfel->lambertianReflectivity;
    const Color3& F_0 = surfel->glossyReflectionCoefficient;
    const Vector3& n  = surfel->shadingNormal;
                    
    float pdfValue = 0.0f;

    // Coarse estimate of Fresnel to influence the sampling pdf only
    const Color3& F = UniversalBSDF::schlickFresnel(F_0, abs(w_o.dot(n)), surfel->smoothness);

    // This only has to be approximate, but must be on [0, 1]
    float probGlossy = F.average() / (F + p_L * (Color3::one() - F)).average();

    // Use a somewhat conservative estimate of the exponent; the true estimator gives values that are too noisy for
    // high smoothness because it doesn't actually match the Disney BRDF. 
//    Vector3::cosPowHemiHemiRandom(w_o.reflectAbout(n), n, UniversalBSDF::smoothnessToBlinnPhongExponent(surfel->smoothness * 0.25f), Random::threadCommon(), w_i, pdfValue);
    Vector3::cosHemiPlusCosPowHemiHemiRandom(w_o.reflectAbout(n), n, UniversalBSDF::smoothnessToBlinnPhongExponent(surfel->smoothness * 0.25f), probGlossy, Random::threadCommon(), w_i, pdfValue);

    const Color3& f = evaluateDisneyBRDF(p_L * (Color3::one() - F_0) + F_0,
        F_0.max(),
        p_L.max(),
        F_0.max(),
        1.0f - surfel->smoothness,
        0.0f, 
        0.0f,
        0.0f,
        0.0f,
        0.0f, 
        0.0f,
        w_i,
        w_o,
        n,
        Vector3::unitX(),
        Vector3::unitZ());

    if (pdfValue < 0.00001f) {
        // Avoid division by zero
        weight = Color3::zero();
    } else {
    	weight =  f * (max(0.0f, w_i.dot(n)) / pdfValue);
    }
}
/////////////////////////////////////////////////////////
// Pete's cone brdf

Color3 evaluateBRDFPeteCone(const Vector3& n, const Vector3& w_o, const Vector3& w_i, const float cos_theta_max) {
    if ((w_i.dot(n) < 0.001f) || (w_o.dot(n) < 0.001f)) {
        // Technically redundant; the cone test includes this case.
        // However, I list it explicitly to match Pete's code
        return Color3::zero();
    }

    const float solid_angle = 2.0f * pif() * (1.0f - cos_theta_max);

    const Vector3& w_mirror = 2.0f * n * n.dot(w_o) - w_o;
    if (w_mirror.dot(w_i) >= cos_theta_max) {
        // In the cone
        return Color3::one() * 1.07f / (solid_angle * max(w_i.dot(n), w_o.dot(n), 0.0f));
        // return Color3::one() / (solid_angle * sqrt(w_i.dot(n), w_o.dot(n), 0.0f));
    } else {
        // Outside of the cone
        return Color3::zero();
    }
}


// Glossy BRDF
void scatterPeteCone(const Vector3& w_o, const Vector3& n, Random& rng, Vector3& w_i, Color3& weight) {
    // Cone angle. This constant is close to mirror reflection:
    const float cos_theta_max = 0.999f;

    const Vector3& w_mirror = 2.0f * n * n.dot(w_o) - w_o;

    // Sample a true cone PDF. For a narrow cone, the PDF will be quite large (100's)
    float pdf_value = 0.0f;

    const bool sampleCone = true;
    if (sampleCone) {
        Vector3::sphericalCapHemiRandom(w_mirror, n, cos_theta_max,rng, w_i, pdf_value);
    } else {
        Vector3::hemiRandom(n, rng, w_i, pdf_value);
    }

    // We didn't importance sample a cosine, so we have to explicitly scale by it here.
    // Divide by the pdf value for the usual Monte Carlo scaling. Because of the cone
    // sampling, the weights should be close to 1.0.
    weight = evaluateBRDFPeteCone(n, w_o, w_i, cos_theta_max) * abs(n.dot(w_i)) / pdf_value;
}

/////////////////////////////////////////////////////////

// Normalized, reciprocal blinn phong (subject to lambertian cosine falloff):  |wh . n|^m * (8 + m) / (8pi)
void scatterBlinnPhong(const Vector3& w_o, const Vector3& n, Random& rng, Vector3& w_i, Color3& weight) {
    float pdf_value;
    Vector3::cosHemiRandom(n, rng, w_i, pdf_value);
    
    // Blinn-Phong Exponent ("shininess")
    const float m = 1000.0f;

    const Color3& f = Color3::one() * pow(max((w_i + w_o).direction().dot(n), 0.0f), m) * ((m + 2.0f) * (m + 4.0f) / (8.0f * pif() * (exp2f(-m / 2.0f) + m)));

    weight = f * abs(n.dot(w_i)) / pdf_value;
}

// Normalized, reciprocal blinn phong with the cone BRDF denominator (resistant to lambertian cosine falloff):  |wh . n|^m * (1 / sqrt(|wi .n * wo . n|)) * (8 + m) / (8pi)
void scatterHackedBlinnPhong(const shared_ptr<UniversalSurfel>& surfel, const Vector3& w_o, Random& rng, Vector3& w_i, Color3& weight) {
    const Vector3& n = surfel->shadingNormal;

    float pdfValue;
    Vector3::cosHemiRandom(n, Random::threadCommon(), w_i, pdfValue);

    // Blinn-Phong Exponent ("shininess")
    const float m = UniversalBSDF::smoothnessToBlinnPhongExponent(surfel->smoothness);

    // supposedly exact normalization: [http://www.thetenthplanet.de/archives/255]
    // ((m + 2.0f) * (m + 4.0f) / (8.0f * pif() * (exp2f(-m / 2.0f) + m)));

    const Color3& f = Color3::one() * pow(max((w_i + w_o).direction().dot(n), 0.0f), m) * ((m + 2.0f) * (m + 4.0f) / (8.0f * pif() * (exp2f(-m / 2.0f) + m))) / sqrt(abs(w_i.dot(n) * w_o.dot(n)));

    weight = f * abs(n.dot(w_i)) / pdfValue;
}

///////////////////////////////////////////////////////////////////////////////////

// http://www.cs.utah.edu/~premoze/dbrdf/dBRDF.pdf
Color3 evaluateDBRDF(const Color3& F0, float smoothness, const Vector3& w_o, const Vector3& n, Vector3& w_i) {
    const Vector3& w_h = (w_o + w_i).direction();
    const Color3& F = UniversalBSDF::schlickFresnel(F0, max(0.0f, w_h.dot(w_i)), smoothness);
    const float m = UniversalBSDF::smoothnessToBlinnPhongExponent(smoothness);
    
    const float nw_i = max(0.0f, w_i.dot(n));
    const float nw_o = max(0.0f, w_o.dot(n));

    // Approximate Phong PDF [not BSDF, thanks to normalization]
    const float p = pow(max(0.0f, n.dot(w_h)), m) * (m + 8.0f);

    return F * (p / (nw_i + nw_o - nw_o * nw_i));
}


void scatterDBRDF(const shared_ptr<UniversalSurfel>& surfel, const Vector3& w_o, Random& rng, Vector3& w_i, Color3& weight) {
	const Vector3& n = surfel->shadingNormal;
	w_i = Vector3::cosHemiRandom(n, rng);

    Color3 f;
    float positiveHemisphere = max(0.0f, sign(w_i.dot(n))) * max(0.0f, sign(w_o.dot(n)));
    const Color3& p_L = surfel->lambertianReflectivity;
    const Color3& F_0 = surfel->glossyReflectionCoefficient;
    weight = evaluateDBRDF(F_0, surfel->smoothness, w_o, n, w_i);
}

///////////////////////////////////////////////////////////////////////////////////



/* Some old GGX code:
	Random& rng = Random::threadCommon();
	const shared_ptr<UniversalSurfel>& us = dynamic_pointer_cast<UniversalSurfel>(surfel);
	const Vector3& n = us->shadingNormal;

	w_i = Vector3::cosHemiRandom(us->shadingNormal, rng);

    Color3 f;
    float positiveHemisphere = max(0.0f, sign(w_i.dot(n))) * max(0.0f, sign(w_o.dot(n)));
    const Color3& p_L = us->lambertianReflectivity;
    const Color3& F_0 = us->glossyReflectionCoefficient;

        ///////////////////////////////////////////////////////
        // GGX
        // From http://graphicrants.blogspot.com/2013/08/specular-brdf-reference.html

        // "Roughness", following UE4 model
        float a = square(1.0f - us->smoothness);
        float a2 = a*a;
        float NoV = max(0.0f, n.dot(w_o));
        float NoL = max(0.0f, n.dot(w_i));
        float G_V = NoV + sqrt( (NoV - NoV * a2) * NoV + a2 );
        float G_L = NoL + sqrt( (NoL - NoL * a2) * NoL + a2 );
        Color3 f_G = F_0 / ( G_V * G_L );


        // Lambertian term
        const Color3 f_L = positiveHemisphere * p_L / pif();

        f = f_L + f_G;
    	weight = pif() * f;
        */