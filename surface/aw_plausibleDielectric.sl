#include <stdrsl/ShadingContext.h>
#include <stdrsl/RadianceSample.h>
#include <stdrsl/SpecularAS.h>
#include <stdrsl/SampleMgr.h>
#include <stdrsl/Lambert.h>
#include <stdrsl/Colors.h>
#include <stdrsl/SphericalHarmonics.h>
#include "../include/aw/aw_TexturingUtils.h"
#include "../include/aw/aw_PrmanEnvLight.h"
/* $Revision: #1 $ $Date: 2011/04/27 $
*/
//  plausibleDielectric:
//  A plausible anisotropic dielectric material. Implements a basic amount 
//  of functionality needed to shade plausible dielectric materials.
//  Uses a stdrsl specular component (stdrsl_SpecularAS) with
//  a simple Lambertian diffuse term, weighted by the transmission coefficient.
class aw_plausibleDielectric(
    uniform float roughness = 0.001;
    uniform float anisotropy = 0;
    uniform float specularGain = 1;
    uniform float diffuseGain = .5;
    uniform color specularColor = color(1);
    uniform color surfaceColor = color(1,1,1);
    uniform float ior = 1.33;
    uniform float mediaIor = 1;
    uniform float displacementAmount = 0;
    uniform string specularMap = "";
    uniform float textureBlur = 1;
    uniform string dispMode = "displace";
    uniform float houdini = 1;
    uniform string surfaceMap = "";
    uniform string roughnessMap = "";
    uniform string displacementMap = "";
    uniform float minSpecularSamples = 1;
    uniform float maxSpecularSamples = 16;
    uniform float indirectDiffuseSamples = 64;
    uniform float indirectDiffuseMaxVar = .02;
    uniform float indirectDiffuseMaxDistance = 1e10;
    uniform string envMap = "";
    uniform string sphericalHarmonicPtc = "";
    uniform string envSpace = "";
    uniform float specularMaxDistance = 1e10;
    uniform float specularBroadening = 1;
    uniform string indirectSpecularSubset = "";
    uniform string lightCategory = "-environment";
    uniform float applySRGB = 1;
    uniform float __computesOpacity = 0;
    uniform float __faceindex = 0;
    output varying color o_indirect_specular = 0;
    output varying color o_diffusecolor = 0;
    output varying color o_specular = 0;
    output varying color o_indirectdiffuse = 0;
    output varying color o_refl_occ = 0;
    output varying color o_occlusion = 0;
    output varying color o_subsurface = 0;
    )
{
    //Member variables
    shader m_lights[];
    shader ptc_lights[];
    stdrsl_ShadingContext m_shadingCtx;
    stdrsl_Fresnel m_fres;
    aw_TexturingUtils m_tex;
    aw_PrmanEnvLight m_pr;
    stdrsl_Lambert m_diffuse;
    stdrsl_SpecularAS m_specular;
    uniform float m_indirectSamples;

    public void construct() {
    }
    public void begin() 
    {
    }
    public void displacement(output point P; output normal N) 
    {
        normal Nn  = normalize(N);
        normal Ngn = normalize(Ng);
        normal normalDelta = Nn - Ngn;
        point Po = transform("object", P);
        m_tex->aw_initFr(textureBlur, houdini, __faceindex);
        if(displacementMap != "" && displacementAmount != 0)
        {
            if(match("ptx$", displacementMap)) {
                color tmp = ptexture(displacementMap, 0, m_tex->m__faceindex, 
                                     m_tex->m_ptxfr);
                m_tex->aw_ptexDisplace(normalize(N), displacementAmount, tmp, dispMode);
            } else {
                color tmp = texture(displacementMap, m_tex->m_fr);
                m_tex->aw_displace(normalize(N), displacementAmount, tmp, dispMode);
            }
        }
        m_shadingCtx->init();
    }
    public void prelighting(output color Ci,Oi) 
    {
        m_lights = getlights("category",lightCategory);
        ptc_lights = getlights("category", "environment");
    }
    //Shared initialization function for diffuse.
    public void initDiff() 
    {
        color surfColor = surfaceColor * Cs;
        if(surfaceMap!= "")
        {
            color c = color texture(surfaceMap, m_tex->m_fr);
            if(applySRGB != 0)
                c = srgbToLinear(c);
            surfColor *= c; 
        }
        m_diffuse->init(surfColor*diffuseGain);
        //Reduce the number of indirect samples based on diffuse ray depth
        stdrsl_SampleMgr mgr;
        mgr->computeDepthBasedSampleReduction(m_shadingCtx->m_DiffuseDepth,
                                              indirectDiffuseSamples, 1,0,
                                              m_indirectSamples);
    }
    //Shared initialization function for specular.
    public void initSpec() 
    {
        m_fres->init(m_shadingCtx, mediaIor, ior);
        color specColor = specularColor; 
        if(specularMap != "")
        {
            color c = color texture(specularMap, m_tex->m_fr);
            if(applySRGB != 0) c = srgbToLinear(c);
            specColor *= c;
        }
        varying float roughscale = 1; if(roughnessMap!= "") { roughscale = 
        texture(roughnessMap[0], m_tex->m_fr); }
        //XXX problems arize when roughness approaches 0
        //roughscale *= m_fres->m_Kt; //reduce roughness at grazing angle
        m_specular->init(m_shadingCtx, specularGain*specColor*m_fres->m_Kr,
                         roughness, anisotropy, roughscale,
                         minSpecularSamples, maxSpecularSamples);
    }
    //Shared initialization function for env light
    public void initEnv()
    {
        color env = 1, ind = 0, indirectlight = 0;
        uniform stdrsl_SphericalHarmonic envmapsh;
        uniform color envmapshcoeffs[];
        stdrsl_SphericalHarmonic dirvissh, radvissh;
        float dirvisshcoeffs[];
        float incradshcoeffs[];
        float ok_occ, ok_rad, exists;
        uniform float ps = arraylength(ptc_lights),i;
        vector m_R = m_fres->m_Rv;
        color albedo = m_diffuse->m_diffColor;
        // if we have aw_prman_envlight
        if(ps >= 1 && textureinfo(sphericalHarmonicPtc, "exists", exists) == 0) {
            for(i = 0; i < ps; i += 1) {
                m_pr->refl_result = ptc_lights[i]->aw_refl_occ(m_R);
                m_pr->occ_result  = ptc_lights[i]->aw_occlusion();
                m_pr->ind_result  = ptc_lights[i]->aw_indirectdiff();
                m_pr->sss_result  = ptc_lights[i]->aw_sss(P);
                
                if(m_pr->refl_result != color(0))
                    m_pr->res += m_pr->refl_result;
                else if(m_pr->occ_result != color(0))
                    m_pr->res +=  m_pr->occ_result;
                if(m_pr->ind_result != color(0))
                    m_pr->res +=  m_pr->ind_result;
                else if(m_pr->sss_result != color(0))
                    m_pr->res += m_pr->sss_result;

                m_pr->indirect = m_pr->ind_result;
            }
        // if we're doing spherical harmonics
        } else if (textureinfo(sphericalHarmonicPtc, "exists", exists) == 1 
                && envMap != "" ) {
            if(textureinfo(envMap, "shcoeffs", envmapshcoeffs) != 0) {
                envmapsh->createFromArray(envmapshcoeffs);
                if (envSpace != "")
                    envmapsh->rotate(envSpace);

                ok_occ = texture3d(sphericalHarmonicPtc, P, m_shadingCtx->m_Ns, 
                        "_dirvisshcoeffs", dirvisshcoeffs);
                ok_rad = texture3d(sphericalHarmonicPtc, P, m_shadingCtx->m_Ns, 
                        "_incradshcoeffs", incradshcoeffs);

                if (ok_occ != 0) {
                    dirvissh->createFromArray(dirvisshcoeffs);
                    m_pr->sph_occ = envmapsh->convolve(dirvissh, 0);
                    m_pr->res += m_pr->sph_occ;
                    m_pr->occ_result = m_pr->sph_occ;
                }
                if (ok_rad != 0) {
                    radvissh->createFromArray(incradshcoeffs);
                    m_pr->sph_rad = envmapsh->convolve(radvissh, 0);
                    m_pr->res += m_pr->sph_rad;
                    m_pr->indirect = m_pr->sph_rad;
                }
            }
        } else { // we trace
            if(indirectDiffuseSamples > 0)  {
                // indirectdiffuse doesn't yet know about m_material
                ind = indirectdiffuse(P, m_shadingCtx->m_Ns,
                                        m_indirectSamples, "adaptive", 1,
                                        "maxvariation", indirectDiffuseMaxVar,
                                        "environmentmap", envMap,
                                        "environmentspace", envSpace,
                                        "environmentcolor", env,
                                        "maxdist", indirectDiffuseMaxDistance);
                if(envMap != "") {
                    m_pr->res += env;
                    m_pr->indirect = env;
                } else {
                    m_pr->res += ind;
                    m_pr->indirect = ind;
                }
            }
        }

    }

    public void lighting(output color Ci, Oi) 
    {
        initDiff();
        initSpec();
        __radiancesample samps[];
        color diff, spec;
        uniform float MISMode = (roughness < .0005) ? 2 : 1; 
        //MISMode = 1;
        directlighting(this, m_lights, "specularresult", spec, 
                        "diffuseresult", diff, "mis", MISMode);
        color indspec = indirectspecular(this, "samplebase", 1,
                                         "subset", indirectSpecularSubset,
                                         "maxdist", specularMaxDistance);
        m_pr->diff_result = diff;
        m_pr->indspec_result = indspec; 
        m_pr->spec_result = spec;

        // prman envlight
        initEnv();
        m_pr->sum();
        o_refl_occ = m_pr->reflAOV;
        o_occlusion = m_pr->occAOV;
        o_subsurface = m_pr->sssAOV;
        o_indirectdiffuse = m_pr->indirect;
        o_specular = m_pr->specAOV;
        o_diffusecolor = m_pr->diffAOV;
        o_indirect_specular = m_pr->indspecAOV;

        diff += m_diffuse->m_diffColor * m_pr->res;
        Ci = diff * m_fres->m_Kt + spec + indspec; // fresnel already in spec
    }
    public void diffuselighting(output color Ci, Oi) 
    {
        // produce view-independent diffuse response (no fresnel)
        initDiff();
        initEnv();
        Ci += directlighting(this, m_lights);
        Ci += m_diffuse->m_diffColor * m_pr->res;
    }
    public void specularlighting(output color Ci,Oi) 
    {
        initSpec();
        color diff = Ci;
        __radiancesample samps[];
        uniform float MISMode = (roughness < .0005) ? 2 : 1; 
        color spec;
        if(specularBroadening == 1)
        {
            varying float tmpRU = m_specular->m_roughValueU;
            varying float tmpRV = m_specular->m_roughValueV;
            m_specular->m_roughValueU = 100;
            m_specular->m_roughValueV = 100;
            spec = directlighting(this, m_lights, "mis", MISMode);
            m_specular->m_roughValueU = tmpRU;
            m_specular->m_roughValueV = tmpRV;
        }
        else
            spec = directlighting(this, m_lights, "mis", MISMode);
        color indspec = indirectspecular(this, "subset", indirectSpecularSubset,
                                         "maxdist", specularMaxDistance);
        Ci = diff * m_fres->m_Kt + spec + indspec;
    }
    public void evaluateSamples(string distribution;
                                output __radiancesample samples[] ) 
    {
        if (distribution == "diffuse")
            m_diffuse->evalDiffuseSamps(m_shadingCtx, m_fres, samples);
        else
            m_specular->evalSpecularSamps(m_shadingCtx, m_fres, samples);
    }
    public void generateSamples(string type; output __radiancesample samples[])
    {
        m_specular->genSpecularSamps(m_shadingCtx, m_fres, type, samples);
    }
}
