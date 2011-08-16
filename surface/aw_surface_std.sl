/************************************************************************
 * aw_surface_std.sl - This is my main uber shader responsible for 90% of
 * the shading in my reel. It's designed to work with aw_prmanEnvlight
 * It's interface is fully represented in a Houdini Digitial Asset.
 *
 * Author: Alan Warren (bluemoonshine@gmail.com)
 *
 * $Revision: 1.1 $    $Date: 2010/02/10 $
 *
 ************************************************************************/


#include "rfm_blinn.h"
#include "normals.h"
#include <stdrsl/SphericalHarmonics.h>
#include <stdrsl/ShadingContext.h>
#include <stdrsl/Fresnel.h>
#include <stdrsl/Math.h>
#include <stdrsl/SampleMgr.h>
#include <stdrsl/Colors.h>
#include "aw/aw_utils.h"
#include "aw/aw_TexturingUtils.h"
#include "aw/aw_PrmanEnvLight.h"

class aw_surface_std(
         uniform string lightcategory = "-environment";
         uniform color surface_id = 1;
         uniform float fresnel_diffuse = 0;
		 /*Koeffs*/
		 uniform float Kd = 0.8, Ks = 0.2, Kr = 0.0, Km = 0, fresnel_reflect = 1;
         uniform color frnl_blend = 1;
         uniform float specmapweight = 0, reflmapweight = 0, texmapweight = 0;
         uniform float diffuseroughness = 0.1;
		 /*texture controls*/
         uniform float maketexlin = 0;
		 uniform string surfacemap = "";
         uniform float blurScale = 1, imagetransparency = 0, alphacomp = 0;
         uniform color background = color(0,0,0), diffusecolor = color(1, 1, 1);
         /* bump map */
         uniform string displacementmap = "";
         uniform string displacemode = "displace";
         /* reflection */
         uniform float makerefllin = 0;
         uniform float refltint = 0;
         uniform float rthreshold = 0.0001;
         uniform float refl_saturation = 0.3;
         uniform float refl_blend = 0.5;
         uniform string reflmap = "";
         uniform float fadereflections = 0;
         uniform float fadescale = 1, fadeexp = 1; // fade ctrl
         uniform float blurEnv = 1;
         uniform float importance_sampling = 0;
         uniform string env_ctx = "all"; // all | surface | env
         /* end switches */
         uniform float samples = 1, angle = 0;
         uniform string envlightcategory = "environment";
         uniform string indenvmap = "";
         uniform float indirectsamples = 0;
         uniform float indirectmaxvar = 0.02;
         uniform float indirectmaxdist = 1e10;
         uniform string envmap = "", envmapspace = "current";
         uniform float yup_to_zup = 0;
		 /*specular*/
         uniform float blinnStrength = 1;
         uniform float is_metallic = 0;
         uniform color hilightcolor = 1;
         uniform float eccentricity = 0.1, rolloff = 0.5;
		 uniform string specmap = "";
         uniform float makespeclin = 0;
		 uniform float nonenv = 0; 
         /* bake to ptex */
         uniform string bakeptexmap = "";
         uniform string ptexmap = "", shptcmap = "";
         uniform float houdini = 1;
         varying float __faceindex = 0;
         // AOV's
		 output varying color o_diffusedirect;
         output varying color o_diffusecolor;
         output varying color o_diffuseindirect;
         output varying color o_refl_occ;
         output varying color o_occlusion;
         output varying color o_subsurface;
         output varying color o_normals;
         output varying color o_indirectdiffuse;
		 output varying color o_specular;
         output varying color o_specularshadow;
         output varying color o_reflection;
         output varying color o_fresnel_kt;
         output varying color o_fresnel_kr;
         output varying color o_id;
         output varying color o_shadow;
         output varying color o_z;
         output varying color o_motion;
         output varying color o_atten;)
{
        shader m_lights[];
        shader env_lights[];
        stdrsl_ShadingContext m_ctx;
        stdrsl_Fresnel m_fres;
        aw_TexturingUtils m_tex;
        aw_PrmanEnvLight m_pr;
        uniform float m_indirectSamples;

        varying normal Nn = 0;
		
		public void construct() {
		}

        public void begin() {
            m_tex->aw_initFr(blurScale, houdini, __faceindex);
        }

        public void prelighting(output color Ci,Oi)
        {
            m_fres->initSchlick(m_ctx, 1, 1.5/*0 ior: no refraction*/, .07);
            env_lights = getlights("category", envlightcategory);
            m_lights = getlights("category", lightcategory);
        }

        public void displacement(output point P; output normal N) 
        {
            normal Nn  = normalize(N);
            normal Ngn = normalize(Ng);
            normal normalDelta = Nn - Ngn;
            point Po = transform("object", P);
            m_tex->aw_initFr(blurScale, houdini, __faceindex);
            if(displacementmap != "" && Km != 0)
            {
                if(match("ptx$", displacementmap)) {
                    color tmp = ptexture(displacementmap, 0, m_tex->m__faceindex, 
                                         m_tex->m_ptxfr);
                    m_tex->aw_ptexDisplace(normalize(N), Km, tmp, displacemode);
                } else {
                    color tmp = texture(displacementmap, m_tex->m_fr);
                    m_tex->aw_displace(normalize(N), Km, tmp, displacemode);
                }
            }
            m_ctx->init();
        }


        public void EnvReflFrame(vector In; normal Nn; float coneSize; float flip; string space;
                              output vector Rn, T0, T1)
        {
           Rn = In - 2 * (In.Nn)*Nn;
           //Rn = reflect(In, Nn);
           Rn = vtransform(space, Rn);
           if(flip != 0) {
               Rn = yup_to_zup(Rn);
           }
           T1 = normalize(Rn ^ Nn);
           T0 = normalize(Rn ^ T1);
           T0 *= coneSize;
           T1 *= coneSize * Nn.Rn;
        }

        public void ReflFrame(vector In; normal Nn; float coneSize; float flip; string space;
                              output vector Rn, T0, T1)
        {
           Rn = In - 2 * (In.Nn)*Nn;

           //Rn = reflect(In, Nn);
           T1 = normalize(Rn ^ Nn);
           T0 = normalize(Rn ^ T1);
           T0 *= coneSize;
           T1 *= coneSize * Nn.Rn;
        }
        private color GatherRefl(point P; vector Rn; vector T0; vector T1; string envmap;
                float flip; string space; uniform float samples; float fade; float fadescale;
                float fadeexp; float env; float importance;, color fresnelmult)
        {
           color ci = 0;
           color hitci = 0;
           color ret = 0;
           vector rdir = 0;
           float raylength = 0;
           // importance samps
           float maxImportance;
           rayinfo("importance", maxImportance);
           float reflSam = Kr * m_fres->m_Kr * samples * maxImportance;
           // normal samps
           uniform float numSamps = 1;
           if(m_ctx->m_RayDepth < 1) {
               numSamps = samples;
           }
           filterregion fr;
           fr->calculate3d(Rn, T0, T1);

           if (importance != 0) {
               if (reflSam > rthreshold) {   // shoot 1 or more rays
                   color weight = color(Kr * m_fres->m_Kr);
                   float sam = max(round(reflSam), 1);
                   ci = 0; hitci = 0;
                   gather("illuminance", P, fr, sam, "weight", weight, 
                          "ray:length", raylength,
                          "surface:Ci", hitci,
                          "ray:direction", rdir) {
                       float fd = fadescale*exp(-fadeexp*raylength);
                       if(fade != 0) {
                            hitci *= fd;
                       }
                       ci += hitci;
                   } else {
                       if(env != 1) {
                          rdir = vtransform(space, rdir);
                          if(flip != 0) 
                              rdir = yup_to_zup(rdir);
                          filterregion frEnv;
                          frEnv->calculate3d(normalize(rdir));
                          frEnv->blur(blurEnv);
                          ci += environment(envmap, frEnv);
                       }
                   }
                    ret += fresnelmult * ci / sam; // mult by fres
               } else if (reflSam / rthreshold > random()) {
                   ci = 0; hitci = 0;
                   gather("illuminance", P, fr, 1,"surface:Ci", hitci) {
                       ci += hitci;
                   } 
                   ret += ci; // no fres mult here to keep weight
               }
               return ret;

           } else {
               ci = 0; hitci = 0;
               gather("illuminance", P, fr, numSamps, "ray:length", raylength, "surface:Ci", hitci,
                       "ray:direction", rdir) {
                   float fd = fadescale*exp(-fadeexp*raylength);
                   if(fade != 0) {
                       hitci *= fd;
                   }
                   ci += hitci;
               } else {
                   if(env != 1) {
                       rdir = vtransform(space, rdir);
                       if(flip != 0) 
                           rdir = yup_to_zup(rdir);
                       filterregion frEnv;
                       frEnv->calculate3d(normalize(rdir));
                       frEnv->blur(blurEnv);
                       ci += environment(envmap, frEnv);
                   }
               }
               return fresnelmult * ci / numSamps;
           }
        }

        private color EnvRefl(vector Rn; vector T0; vector T1; string envmap;)
        {
            color ci = 0;
            filterregion frEnv;
            frEnv->calculate3d(Rn, T0, T1);
            frEnv->blur(blurEnv);

	        if(envmap != "") { 
                ci = color environment(envmap, frEnv);
            }
            return ci;
        }

        private color FrnlBlend(color vNormal; color vGrazing; float cosTheta)
        {  // schlick approx
           float wt = pow(1-cosTheta, 5);
           return mix(vNormal, vGrazing, color(wt));
        }

        public void initEnv()
        {
            color env = 1, ind = 0, indirectlight = 0;
            uniform stdrsl_SphericalHarmonic envmapsh;
            uniform color envmapshcoeffs[];
            stdrsl_SphericalHarmonic dirvissh, radvissh;
            float dirvisshcoeffs[];
            float incradshcoeffs[];
            float ok_occ, ok_rad, exists;
            uniform float ps = arraylength(env_lights),i;
            vector m_R = m_fres->m_Rv;
            color albedo = m_pr->color_result;
            // if we have aw_prman_envlight
            if(ps >= 1 && textureinfo(shptcmap, "exists", exists) == 0) {
                for(i = 0; i < ps; i += 1) {
                    if(env_lights[i] != null) {
                        m_pr->refl_result = env_lights[i]->aw_refl_occ(m_R);
                        m_pr->glossy_result = env_lights[i]->aw_refl_indirect(m_R);
                        m_pr->occ_result  = env_lights[i]->aw_occlusion();
                        m_pr->ind_result  = env_lights[i]->aw_indirectdiff();
                        m_pr->sss_result  = env_lights[i]->aw_sss(P);
                        m_pr->res = color(1);
                        
                        if(m_pr->refl_result != color(0))
                            m_pr->res *= m_pr->refl_result;
                        //else if(m_pr->glossy_result != color(0))
                            //m_pr->res *= m_pr->glossy_result;
                        else if(m_pr->occ_result != color(0))
                            m_pr->res *=  m_pr->occ_result;
                        if(m_pr->ind_result != color(0))
                            m_pr->indirect +=  m_pr->ind_result;
                        else if(m_pr->sss_result != color(0))
                            m_pr->sss += m_pr->sss_result;
                    }
                }
            // if we're doing spherical harmonics
            } else if (textureinfo(shptcmap, "exists", exists) == 1 
                    && indenvmap != "" ) {
                if(textureinfo(indenvmap, "shcoeffs", envmapshcoeffs) != 0) {
                    envmapsh->createFromArray(envmapshcoeffs);
                    if (envmapspace != "")
                        envmapsh->rotate(envmapspace);

                    ok_occ = texture3d(shptcmap, P, m_ctx->m_Ns, 
                            "_dirvisshcoeffs", dirvisshcoeffs);
                    ok_rad = texture3d(shptcmap, P, m_ctx->m_Ns, 
                            "_incradshcoeffs", incradshcoeffs);
                    m_pr->res = color(1);
                    if (ok_occ != 0) {
                        dirvissh->createFromArray(dirvisshcoeffs);
                        m_pr->sph_occ = envmapsh->convolve(dirvissh, 0);
                        m_pr->res *= m_pr->sph_occ;
                        m_pr->occ_result = m_pr->sph_occ;
                    }
                    if (ok_rad != 0) {
                        radvissh->createFromArray(incradshcoeffs);
                        m_pr->sph_rad = envmapsh->convolve(radvissh, 0);
                        m_pr->res *= m_pr->sph_rad;
                        m_pr->indirect = m_pr->sph_rad;
                    }
                }
            } else { // we trace
                if(indirectsamples > 0)  {
                    m_pr->res = color(1);
                    ind = indirectdiffuse(P, m_ctx->m_Ns,
                                            m_indirectSamples, "adaptive", 1,
                                            "maxvariation", indirectmaxvar,
                                            "environmentmap", indenvmap,
                                            "environmentspace", envmapspace,
                                            "environmentcolor", env,
                                            "maxdist", indirectmaxdist);
                    if(indenvmap != "") {
                        m_pr->res += env;
                        m_pr->indirect = env;
                    } else {
                        m_pr->res += ind;
                        m_pr->indirect = ind;
                    }
                }
            }

        }

        public void getDiffuse(output color di, diffusedirect, occlusiondirect, specC, specK; float roughness)
            {
                color k = color(0);
                color cldiff = color(0);
                color inshadow, noshadow;
                float C1,C2,C3,ff,fP,L1,L2,sigma,A,B,theta_r;
                color CC =0;
                vector VpN;
        
                normal Ns = shadingnormal(N);
                vector V = normalize(-I);
                float ndotv = Ns.V;
        
                sigma = roughness * roughness; 
                A = 1 - .5 * sigma / (sigma + .33);
                B = .45 * sigma / (sigma + .09);
                theta_r = acos(ndotv);
                VpN = normalize(V-Ns*(ndotv));
                //specular init
                color highlights_component = 0;
                color hcolor = (is_metallic == 1) ?
                    m_pr->color_result*hilightcolor : hilightcolor; 

                float computedEcc = eccentricity * eccentricity - 1;
                if(computedEcc > -.00001 && computedEcc < .00001)
                    computedEcc = .00001;
                computedEcc = 1.0 / computedEcc;
                //diffuse
                 sigma = roughness * roughness; 
                 A = 1 - .5 * sigma / (sigma + .33);
                 B = .45 * sigma / (sigma + .09);
                 theta_r = acos(ndotv);
                 VpN = normalize(V-Ns*(ndotv));

                 // diffuse
                 illuminance("-environment", P, Ns, 1.57)
                 {
                    vector Ln = normalize(L);
                    float cos_theta_i = Ln.Ns;
                    float cos_phi_d = VpN.(Ln - Ns*cos_theta_i);
		            float theta_i = acos(cos_theta_i);
		            float alpha = max(theta_i, theta_r);
		            float beta = min (theta_i, theta_r);
		            C1 = A;
		            if (cos_phi_d >= 0) {
		                 C2 = B;
		            } else {
		                 ff = (2*beta)/ PI;
		                 fP = ff*ff;
		                 C2 = .45 *(sigma/(sigma+.09)) * (sin(alpha)- (pow(ff,3)));
                    }
                    C3 = 0.125 * (sigma/(sigma+0.09)) * pow((4*beta*alpha)/(PI*PI),2);
                    L1 = (cos_theta_i * (C1 + cos_phi_d * C2 * tan(beta) + (1- abs(cos_phi_d)) * 
                          C3 * tan((alpha+beta)/2)));
                    L2 = 0.17*(cos_theta_i)*(sigma/(sigma+0.13)) * (1 - (cos_phi_d *pow((2*beta/PI),2) ));

		         	
		            float nondiff = 0;
		            lightsource("__nondiffuse", nondiff);
		            if (nondiff < 1)
		            {
		                float k = (1-nondiff) * (L1 + L2);
		                if( 0 != lightsource("_shadow", inshadow) )
		                     occlusiondirect += inshadow;
		                if( 0 == lightsource("_cl_noshadow", noshadow) )
		                     noshadow = Cl;
		                cldiff = noshadow - Cl;
		                diffusedirect += k * noshadow; //directC
		                di += k * cldiff; //lightC
		            }

                    // specular
                    float nonspec = 0;
                    lightsource("__nonspecular", nonspec);
                    if (nonspec < 1)
                    {
                        rfmBlinnSpecular(Ln, Cl, noshadow, 1-nonspec,
                                m_ctx->m_Ns, m_ctx->m_Vn, ndotv, computedEcc, rolloff,
                                blinnStrength, m_ctx->m_RayDepth,
                                specC, specK);
                        specC *= hcolor;
                    }
		         } // illuminance() 

            }

		public void lighting(output color Ci, Oi) {
			float scos, ssin, dot_i, maxrim, nhits =0;
			vector lightdir;
			color texcol = 1, colcomp = 1, 
			specweight = 1, shiny = 0, hitc = 0, rt = color(0);
            vector raydir = 0;
			vector Ln, H, D;
            color specK=0, C=0;
            float opacity = 1;

            ///////////////////////////////////////////////////////////////////////////
			/* texture map                                                         'b*/
            ///////////////////////////////////////////////////////////////////////////

			if(surfacemap != "") {
				if(match("ptx$", surfacemap)) {
					texcol = ptexture(surfacemap, 0, m_tex->m__faceindex, m_tex->m_ptxfr);
				} else { // normal texture map
        	    	opacity = texture(surfacemap[3], m_tex->m_fr);
					texcol = texture(surfacemap, m_tex->m_fr);
        	    	texcol = mix(texcol, color(1), texmapweight);
        	    	if(alphacomp != 0) {
        	    		texcol = mix(texcol, background, 1-opacity);
        	    	}
				}
			}

            if(maketexlin != 0) 
                texcol = srgbToLinear(texcol);
            //push color to diffuse struct

            //Oi = Os;
            //if(imagetransparency != 0) {
            //   Oi = Os * mix(color(1), color(0), 1-opacity);
            //}

            // sum colors
            color surfcolor = Kd * (Cs * diffusecolor * texcol);
            m_pr->color_result = surfcolor;

            ///////////////////////////////////////////////////////////////////////////
            /* shadow pass */
            ///////////////////////////////////////////////////////////////////////////
            illuminance("shad", P) {
            		varying color _shadow;
                    varying color shad_mult;
                    varying float _atten;
                    varying color o_atten = 1;
            		uniform float li_disableshadpass;
            	    // import shadow	
            		lightsource("_shadow", _shadow);
                    // bring in attenuation without any decay
                    lightsource("_atten", _atten);
            		lightsource("disableshadpass", li_disableshadpass);
    
            		if(li_disableshadpass != 1){
            			shad_mult = (_shadow[0]+_shadow[1]+_shadow[2])/3;
                        shad_mult = mix(color(0), color(1), shad_mult);
                        //take attenuation and clamp it
                        _atten = clamp(_atten,0,1);
                        o_atten = color(_atten,_atten,_atten);
                        //incorporate attenuation with shadow
                        //shad_mult *= o_atten;
            		} else {
            				shad_mult = 1;
            		}
                    m_pr->shadow_result = shad_mult;
            }


            ///////////////////////////////////////////////////////////////////////////
			/* diffuse + specular illumination                                     'c*/
            ///////////////////////////////////////////////////////////////////////////
            color di = color(0);
            color diffusedirect = color(0);
            color occlusiondirect = color(0);
            getDiffuse(di, diffusedirect, occlusiondirect, C, specK, diffuseroughness );
            m_pr->diffdirect_result = diffusedirect;

            // depth based sample reduction for traced indirectdiffuse only
            stdrsl_SampleMgr mgr;
            mgr->computeDepthBasedSampleReduction(m_ctx->m_DiffuseDepth,
                                                  indirectsamples, 0, 1,
                                                  m_indirectSamples);

            ///////////////////////////////////////
            initEnv(); // initialize environment //
            ///////////////////////////////////////

            ///////////////////////////////////////////////////////////////////////////
            /* reflection                                                          'f*/
            ///////////////////////////////////////////////////////////////////////////
            vector Rn, T0, T1;
            vector En, T2, T3;
	        color Cr = 0;
            color rho = 0;
            vector In = normalize(I);
            normal Nn = normalize(N);

            EnvReflFrame(In, Nn, angle, yup_to_zup, envmapspace, En, T2, T3);
            ReflFrame(In, Nn, angle, yup_to_zup, envmapspace, Rn, T0, T1);
            rho = FrnlBlend(frnl_blend, color(1), abs(m_ctx->m_In.m_ctx->m_Ns));

            if(Kr > 0 && env_ctx == "all" ) {
                rho *= fresnel_reflect;
                rt += GatherRefl(P, Rn, T0, T1, envmap, yup_to_zup, envmapspace, 
                        samples, fadereflections, fadescale, fadeexp, 0, importance_sampling, rho );

            } else if(Kr > 0 && env_ctx == "env")  {
	             rt += EnvRefl(En, T2, T3, envmap);
             
            } else if(Kr > 0 && env_ctx == "surface") {
                rho *= fresnel_reflect;
                rt += GatherRefl(P, Rn, T0, T1, envmap, yup_to_zup, envmapspace, 
                        samples, fadereflections, fadescale, fadeexp, 1, importance_sampling, rho );
            } else if(Kr > 0 && env_ctx == "pointcloud") {
                rt += m_pr->glossy_result;
            }

            color reflhsv = 0;
            color reflrgb = 0;
            color myrt = rt;
            reflhsv = ctransform("hsv", myrt);

            if(refltint != 0) {
                setcomp(reflhsv, 1, refl_saturation);
                reflrgb = ctransform("hsv", "rgb", reflhsv);
                rt = reflrgb;
                rt = mix(rt, surfcolor, refl_blend);
            }

            /* reflection map */ 
            color Rmap = 1;
            if(reflmap != "")
            {
			   if(match("ptx$", reflmap)) {
					Rmap = ptexture(reflmap, 0, m_tex->m__faceindex, m_tex->m_ptxfr);
			   } else {
               	  Rmap = texture(reflmap, m_tex->m_fr);
               	  Rmap = mix(Rmap, color(1), reflmapweight);
               	  // linear map
               	  if(makerefllin != 0)
               	  {
               	     Rmap = srgbToLinear(Rmap);
               	  }
            	}
			}
            m_pr->gather_result = Kr * (Rmap * rt);


            // direct + indirect lighting
            color justdiff, justfres;
            justdiff = (diffusedirect - di) + m_pr->indirect;
            // store diffuse result
            m_pr->diff_result = justdiff;
            m_pr->diffind_result = justdiff;
            // store fresnel calcs
            m_pr->fresKt_result = m_fres->m_Kt*color(1);
            m_pr->fresKr_result = m_fres->m_Kr*color(1);
            // multiply diffuse with surface color
			color illumsurfcolor = m_pr->diff_result * surfcolor;
            // fresnel
            justfres = fresnel_diffuse * m_fres->m_Kt;
            if (fresnel_diffuse > 0)
                illumsurfcolor *= justfres;

            /* spec map */	
			if(specmap != "") {
			   if(match("ptx$", specmap)) {
					specweight = ptexture(specmap, 0,m_tex->m__faceindex, m_tex->m_ptxfr);
			   } else {
			 		specweight = texture(specmap, m_tex->m_fr);
             		specweight = mix(specweight, color(1), specmapweight);

             		if(makespeclin != 0) {
             		       specweight = srgbToLinear(specweight);
             	    }
				} 
			}

			m_pr->spec_result = (Ks*specweight) * C;
            m_pr->specshad_result = (Ks*specweight) * specK;

            // holdout coshaders
             shader matteCo[] = getshaders("category", "holdouts");
             uniform float i,matteN = arraylength(matteCo);

             for (i = 0; i < matteN; i += 1)
             {
                 if (matteCo[i] != null)
                 {
                     matteCo[i]->out();
                 }
             }
				
            // normal output
        	float nmapx = 0.5*(1+xcomp(m_ctx->m_Ns));
        	float nmapy = 0.5*(1+ycomp(m_ctx->m_Ns));
        	float nmapz = 0.5*(1+zcomp(m_ctx->m_Ns));
			float a = area(P, "dicing");
        	m_pr->normal_result = color(nmapx, nmapy, nmapz);

            vector motion = transform("NDC", dPdtime);
            float tmapx = 0.5*(1+xcomp(motion));
            float tmapy = 0.5*(1+ycomp(motion));
            float tmapz = 0.5*(1+zcomp(motion));

            //depth
            float zDepth = 1 - distance(E,P);

            // aov's
            m_pr->sum();
            o_diffusedirect = m_pr->diffdirAOV;
            o_diffusecolor = m_pr->colorAOV;
            o_indirectdiffuse = ((1-nonenv) * m_pr->indirect);
            o_subsurface = ((1-nonenv) * m_pr->sss);
            o_refl_occ = ((1-nonenv) * m_pr->reflAOV);
            o_reflection = m_pr->gathAOV;
            o_occlusion = ((1-nonenv) * m_pr->occAOV);
            o_specular = m_pr->specAOV;
            o_specularshadow = m_pr->specshadAOV;
            o_shadow = m_pr->shadowAOV;
            o_normals = m_pr->normalAOV;
            o_fresnel_kt = m_pr->fresKtAOV;
            o_fresnel_kr = m_pr->fresKrAOV;
            o_diffuseindirect = m_pr->diffindAOV;
            o_motion = color(tmapx, tmapy, tmapz);
            o_id = surface_id;
            o_z = mix(color(0), color(1), zDepth);

            /* Grand Total */
            color finalsum = illumsurfcolor + m_pr->gather_result + m_pr->spec_result;

            ///////////////////////////////////////////////////////////////////////////
            /* bake illuminance to ptex                                            'g*/
            ///////////////////////////////////////////////////////////////////////////

            if(m_ctx->m_RayDepth < 1 && bakeptexmap != "") {

                bake3d(bakeptexmap, "illumination", m_tex->m_facedata, m_ctx->m_Ns, "coordsystem",
                        "_disable", "interpolate", 1,
                        "illumination", finalsum);
            }

            color ptxdata = 0;
            if(ptexmap != "") {
               ptxdata = ptexture(ptexmap, 0, m_tex->m__faceindex, m_tex->m_ptxfr); 
			   finalsum = ptxdata;
            } 

			Oi = 1;
		    Ci = finalsum * Oi;

		}
}

	
