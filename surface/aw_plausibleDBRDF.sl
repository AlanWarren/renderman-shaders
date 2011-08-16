/* $Revision: #1 $ $Date: 2011/04/27 $
# ------------------------------------------------------------------------------
#
# Copyright (c) 2011 Pixar Animation Studios. All rights reserved.
#
# The information in this file is provided for the exclusive use of the
# software licensees of Pixar.  It is UNPUBLISHED PROPRIETARY SOURCE CODE
# of Pixar Animation Studios; the contents of this file may not be disclosed
# to third parties, copied or duplicated in any form, in whole or in part,
# without the prior written permission of Pixar Animation Studios.
# Use of copyright notice is precautionary and does not imply publication.
#
# PIXAR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING
# ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT
# SHALL PIXAR BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES
# OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
# WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
# ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS
# SOFTWARE.
#
# Pixar
# 1200 Park Ave
# Emeryville CA 94608
#
# ------------------------------------------------------------------------------
*/

#include <stdrsl/ShadingContext.h>
#include <stdrsl/RadianceSample.h>
#include <stdrsl/Math.h>
#include <stdrsl/SpecularDBRDF.h>
#include <stdrsl/SampleMgr.h>
#include <stdrsl/Lambert.h>
#include <stdrsl/Colors.h>
#include <stdrsl/Fresnel.h>
#include <stdrsl/SphericalHarmonics.h>
#include "aw/aw_TexturingUtils.h"
#include "aw/aw_PrmanEnvLight.h"

//
// helloDBRDF:
//  the hello world of arealight-aware shaders using a stdrsl specular
//  component (stdrsl_SpecularDBRDF) with a simple Lambertian diffuse term.

class aw_plausibleDBRDF(
        uniform string distribution = "beckmann";
        uniform color diffuseGain = .5;
        uniform color diffuseColor = color(.18);
        uniform float specularGain = 1;
        uniform color specularColor = color(1);
        uniform float specularRoughness = .01;
        uniform float specularSamples = 16;
        uniform float displacementAmount = 0;
        uniform float textureBlur = 1;
        uniform string dispMode = "displace";
        uniform float houdini = 1;
        uniform string surfaceMap = "";
        uniform string specularMap = "";
        uniform string reflectionMap = "";
        uniform string displacementMap = "";
        uniform float indirectDiffuseSamples = 64;
        uniform float indirectDiffuseMaxVar = .02;
        uniform float indirectDiffuseMaxDistance = 1e10;
        uniform string envMap = "";
        uniform string sphericalHarmonicPtc = "";
        uniform string envSpace = "";
        uniform float specularMaxDistance = 1e10;
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
    shader m_lights[];
    shader ptc_lights[];
    stdrsl_ShadingContext m_shadingCtx;    
    stdrsl_SpecularDBRDF m_specular;
    stdrsl_Lambert m_diffuse;
    stdrsl_Fresnel m_fres;
    aw_TexturingUtils m_tex;
    aw_PrmanEnvLight m_pr;
    uniform float m_indirectSamples;

    public void construct()
    {
    }

    public void begin()
    {
        m_shadingCtx->init();
        m_fres->initSchlick(m_shadingCtx, 1, 0, .07);
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
    //Shared initialization function for diffuse.
    public void initDiff() 
    {
        color surfColor = diffuseColor * Cs;
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
    public void initSpec()
    {
        color spec = specularGain * specularColor * m_fres->m_Kr;
        if(specularMap != "")
        {
            color c = color texture(specularMap, m_tex->m_fr);
            if(applySRGB != 0) c = srgbToLinear(c);
            spec *= c;
        }

        m_specular->init(m_shadingCtx,
                         distribution, spec, specularRoughness, 
                         1/*minsamps*/, specularSamples);
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
    public void prelighting(output color Ci, Oi)
    {
        //color spec = specularGain * specularColor * m_fres->m_Kr;
        //m_specular->init(m_shadingCtx,
                         //distribution, spec, specularRoughness, 
                         //1/*minsamps*/, specularSamples);
        m_lights = getlights("category",lightCategory);
        ptc_lights = getlights("category", "environment");

    }

    public void lighting(output color Ci, Oi)
    {
        initDiff();
        initSpec();
        // In lighting we deliver full solution, obtain the diffuse and
        // specular results and apply fresnel to diffuse term.
        __radiancesample samps[];
        color diff, spec;        
        uniform float MISMode = (specularRoughness < .0005) ? 2 : 1;

        directlighting(this, m_lights,
                       "specularresult", spec, 
                        "diffuseresult", diff,
                        "mis", MISMode,
                        "materialsamples", samps);
        color indspec = indirectspecular(this, "materialsamples", samps,
                                         "subset", indirectSpecularSubset,
                                         "maxdist", specularMaxDistance);
        m_pr->diff_result = diff;
        m_pr->spec_result = spec;
        m_pr->indspec_result = indspec;

        // prman_envlight
        initEnv();
        m_pr->sum();
        o_refl_occ = m_pr->reflAOV;
        o_occlusion = m_pr->occAOV;
        o_subsurface = m_pr->sssAOV;
        o_indirectdiffuse = m_pr->indirect;
        o_diffusecolor = m_pr->diffAOV;
        o_specular = m_pr->specAOV;
        o_indirect_specular = m_pr->indspecAOV; 

        diff += m_diffuse->m_diffColor * m_pr->res;
        Ci = diff * m_fres->m_Kt + spec + indspec;
    }

    public void diffuselighting(output color Ci, Oi)
    {
        // produce view-independent diffuse response (no fresnel)
        initDiff();
        initEnv();
        Ci += directlighting(this, m_lights);
        Ci += m_diffuse->m_diffColor * m_pr->res;
    }

    public void specularlighting(output color Ci, Oi)
    {
        // in specularlighting we produce final radiance assuming
        // diffuse result is present in Ci. We share  samples generated
        // by indirectspecular with directlighting.
        initSpec();
        color diff = Ci;
        __radiancesample samps[];
        uniform float MISMode = (specularRoughness < .0005) ? 2 : 1;
        color spec = directlighting(this, m_lights, "mis", MISMode, 
                                    "materialsamples", samps);
        color indspec = indirectspecular(this, "materialsamples", samps,
                                         "subset", indirectSpecularSubset,
                                         "maxdist", specularMaxDistance);
        Ci = diff * m_fres->m_Kt + spec + indspec;
    }


    /* integrator callbacks ----------------------------------------------*/

    public void evaluateSamples(string distribution;
                                output __radiancesample samples[])
    {
        if(distribution == "diffuse")
            m_diffuse->evalDiffuseSamps(m_shadingCtx, m_fres, samples);
        else
        if(distribution == "specular")
            m_specular->evalSpecularSamps(m_shadingCtx, m_fres, samples);
    }

    // generateSamples must fill in the direction and distance (inf) for 
    // each samples, as well as the materialResponse (brdf value) and 
    // the materialPdf
    public void generateSamples(string type;
                                output __radiancesample samples[])
    {
        if(type == "specular")
        {
            m_specular->genSpecularSamps(m_shadingCtx, m_fres, type, samples);
        }
    }
}
