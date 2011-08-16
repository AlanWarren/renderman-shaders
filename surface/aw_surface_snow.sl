/************************************************************************
 * aw_surface_snow.sl - This is my main uber shader, but with snow based
 * on surface normals added. Uses Buratti diffuse for snow portions.
 * note: this has not been converted to RSL 2.0, and is thus a bit slower
 *
 * Author: Alan Warren (bluemoonshine@gmail.com)
 *
 * $Revision: 1.1 $    $Date: 2010/02/10 $
 *
 ************************************************************************/

#include "normals.h"
#include "aw_utils.h"
#include "rfm_blinn.h"
#include "pxslRayUtil.h"
#include "pxslUtil.h"
#include "filterwidth.h"
#include "shrimp/rsl_shrimp_fractal.h"
#include "shrimp/rsl_shrimp_shadingmodels.h"
#include <stdrsl/SphericalHarmonics.h>

class aw_surface_snow (
		 /*Koeffs*/
		 float Ka = 1, Kd = 0.8, Ks = 0.2, Kr = 0.0, Km = 0, Krim = 0, Kb = 0, frnl = 1;
         float smapw = 0, rmapw = 0, tmapw = 0;
		 color snowcolor = color(1,1,1);
         float reflectionSpecularity = 1;
		 string dirtmap = "";
		 float displace_all = 0;
		 string dispspace = "current";
		 float snow_ht = 0.5;
		 float snowK = 0.5, sden = 1;
		 float snow_noise_scale = 0.1;
		 float snow_fresenel_coeff = 0; // 0 - 1
		 float snow_height_ratio = 0.5;
		 vector driftDir = vector(0,1,0);
		 // noise
		 float noise_freq = 12;
		 float disp_freq = 4;
		 float non_snow_disp = 0.2;
		 float dispscale = 0;
		 float pockmarks = 0.1;
		 float noise_amp = 1.2;
		 float fwidth_mult = 0.05;
		 float octaves = 6;
		 float lacunarity = 2;
		 float gain = 0.6;
	     float K_angle_mask = 0.5;	
		 float spark_spec = 1;
		 float spark_roughness = 0.1;
		 float spark_freq = 24;
		 float spark_amp = 1;
		 float spark_octaves = 6;
		 color sparkle1 = color(1,1,1);
		 color sparkle2 = color(1,1,1);
		 float spark_gain = 0.6;
		 float spark_lacunarity = 2;
		 float spark_pow = 1;
		 color add_spec_color = color(0.74, 0.91, 1);
		 float Ks_add = 0.7;
	     float add_spec_roughness = 0.11;
		 /*texture controls*/
         float maketexlin = 0;
		 string mapname = "";
         float srep = 1,trep = 1;
         float blurScale = 1, imagetransparency = 0, alphacomp = 0;
         color background = color(0,0,0), diffusecolor = color(1, 1, 1);
         /* bump map */
         float makebumplin = 0;
         string texbump = "";
		 float bumpnoise = 0;
         /* reflection */
         float makerefllin = 0;
         float refltint = 0;
		 float graze = 4;;
         float refl_saturation = 0.3;
         float refl_blend = 0.5;
         color frnl_blend = 1; 
		 string reflmap = "";
         float fadereflections = 0;
         float fadescale = 1, fadeexp = 1; // fade ctrl
         float blurEnv = 1;
         string env_ctx = "all"; // all | surface | env
         /* end switches */
         float samples = 1, angle = 0;
         string envmap = "", envmapspace = "current";
         float yup_to_zup = 0;
		 /*rim light*/
		 color rimcolor = color(1,1,1), lightcolor = 1;
         float rimstrength = 0, rimthreshold = 0.1;
		 /* sss */
		 string sssmap = "";
		 color albedo = color(0.830, 0.791, 0.753);
		 color dmfp = color(8.51, 5.57, 3.95); // marble
		 float ior = 1.5;
		 float unitlength = 0.1;
		 float smooth = 0;
		 float maxsolidangle = 1.0; // quality knob. lower = better
		 /*specular*/
         float specmode = 0;
         float anisotropy = 1;
         float specular_gloss = 1;
         float is_metallic = 0;

         color hilightcolor = 1;
         float eccentricity = 0.1, rolloff = 0.5;
		 string specmap = "";
         float makespeclin = 0;
		 float nonenv = 0; 
         /* bake to ptex */
         string bakemap = "", bakedata = "occlusion", bakemode = "ptex";
         string displaychannel = "occlusion";
         float interpolate = 1, bakenormal = 1, imagebake = 0;
         string ptexmap = "", ptcmap = "";
         varying float __faceindex = 0;
		 output varying color o_diffusedirect;
         output varying color o_diffusecolor;
         output varying color o_diffuseindirect;
		 output varying color o_spark;
		 output varying color o_add_spec;
         output varying color o_refl_occ;
         output varying color o_occlusion;
         output varying color o_subsurface;
         output varying color o_normals;
         output varying color o_colorbleeding;
		 output varying color o_specular;
         output varying color o_specularshadow;
		 output varying color o_rimcolor;
         output varying color o_reflection;
         output varying color o_shadow;)
{
		constant float lerp;
        ShadingContext m_shadingCtx;

        varying normal Nn = 0;
		varying float hump = 0;
		varying float snowmask = 0;
		uniform string sys = "";


		public void construct() {
			lerp = 1;
		}

        public void begin() {
            m_shadingCtx->InitializeAttributes();
        }

        public void EnvReflFrame(vector In; vector Nn; float coneSize; float flip; string space;
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

        public void ReflFrame(vector In; vector Nn; float coneSize; float flip; string space;
                              output vector Rn, T0, T1)
        {
           Rn = In - 2 * (In.Nn)*Nn;

           //Rn = reflect(In, Nn);
           T1 = normalize(Rn ^ Nn);
           T0 = normalize(Rn ^ T1);
           T0 *= coneSize;
           T1 *= coneSize * Nn.Rn;
        }

        private color GatherRefl(point P; vector Rn; vector T0; vector T1; string envmap; float flip; string space; 
                uniform float samples; float fade; float fadescale; float fadeexp; float env;)
        {
           color ret = 0;
           color Cret = 0;
           vector rdir = 0;
           float raylength = 0;

           uniform float numSamps = 1;
           if(m_shadingCtx->m_raydepth < 1) {
               numSamps = samples;
           }

           filterregion fr;
           fr->calculate3d(Rn, T0, T1);

           gather("illuminance", P, fr, numSamps, "ray:length", raylength, "surface:Ci", Cret, 
                   "ray:direction", rdir) {
               float fd = fadescale*exp(-fadeexp*raylength);
               if(fade != 0) {
                    Cret *= fd;
               }
               ret += Cret;
           } else {
                rdir = vtransform(space, rdir);
                if(flip != 0) rdir = yup_to_zup(rdir);
              if(env != 1) {
                filterregion frEnv;
                frEnv->calculate3d(normalize(rdir));
                //yup_to_zup(rdir);
                frEnv->blur(blurEnv);
                ret += environment(envmap, frEnv);
                //ret += ptexture(envmap, 0, __faceindex, frEnv);
               }
           }
           return ret / numSamps;
        }

        private color EnvRefl(vector Rn; vector T0; vector T1; string envmap;)
        {
            color ret = 0;
            filterregion frEnv;
            frEnv->calculate3d(Rn, T0, T1);
            frEnv->blur(blurEnv);

	        if(envmap != "") { 
                ret = color environment(envmap, frEnv);
            }
            return ret;
        }

        private color FrnlBlend(color vNormal; color vGrazing; float cosTheta)
        {  // schlick approx
           float wt = pow(1-cosTheta, 5);
           return mix(vNormal, vGrazing, color(wt));
        }

        private float roughness(float gloss)
        {
            return pow(2.0, 8.0 * gloss);
        }


        // Compute sampling parameters.
        private float compute_sampling_angle( float glossiness; )
        {
        	if( glossiness >= 1.0 )
        		return 0.0;
        	else
        
        	    return (1.0 - pow(clamp(glossiness, 0.0, 1.0), 0.2)) * PI / 2.0;
        }

		private float spark_noise(point p; float amp; float octaves; float ain; float lacunarity; float fwidth;)
		{
			float sum = 0;
			point pp = p;
			float fw = fwidth;
			uniform float i;

			for (i = 0;  i < octaves;  i += 1) {
			#pragma nolint
				sum += amp * gain * filteredsnoise(p, fw);
				//amp *= gain;
				pp *= lacunarity;
				fw *= lacunarity;
			}	
			
			return sum;
		}

		private color category_specular( normal Nf; vector V; float rolloff;) {
			color C = 0;
			extern color Cl;
			extern point P;
		
			illuminance( "-environment", P, Nf, PI/2 ) {
				float nonspec = 0;
				lightsource ("__nonspecular", nonspec);
		
				C += Cl * (1 - nonspec) * specularbrdf (normalize(L), normalize(Nf), normalize(V), rolloff);
			}
		
			return C;
		}

		private void pxr_brownian(uniform float baseFreq, octaves, lacunarity, gain, rescale;
		                 float fourthdimension;
						 varying point p;
						 vector dpu, dpv;
						 output float result;)
		 {
				 uniform float amp = 1;
				 varying point pp = baseFreq * p;
				 varying float fw = pxslFilterWidth(baseFreq*dpu, baseFreq*dpv);
				 uniform float i;
				 result = 0;
				 
				 for (i = 0; i < octaves; i+=1) {
		#pragma nolint
						 result += amp * filteredsnoiset(pp, fw, fourthdimension);
						 amp *= gain; pp *= lacunarity; fw *= lacunarity;
				 }
				 if (rescale != 0)
				 {
						 result = .5 * (1 + result);
						 if(rescale == 2)
						 {
								 result = clamp(result, 0, 1);
						 }
				 }
		 }

		private void getSnow(varying normal dir; uniform float snowline;  uniform vector snowDrift;
							 output varying float c;)
		{
			vector snowDir = transform("object", normalize(snowDrift));
			vector objectDir = transform("object", normalize(dir));
			float angle = objectDir.snowDir;
			//angle = clamp(angle, 0, 1);
			//angle = pow(angle, snowexponent);

			if(angle >= snowline)
				c = 1.0;
			else
				c = 0.0;
		}

		private void getSnow2(varying normal dir; uniform float snowline; varying float altitude;
							 output varying float c;)
		{
			vector objectDir = transform("object", normalize(dir));

			if(altitude >= snowline)
				c = 1.0;
			else
				c = 0.0;

			c *= 1 - smoothstep(ycomp(objectDir)-0.1,ycomp(objectDir)+0.1,snow_height_ratio);
		}

		public void displacement(output point P; output  normal N) {
			Nn = normalize(N);
			vector In = normalize(I);
			normal NG = normalize(Ng);

			if(dispspace == "")
				sys = "shader";
			else
				sys = dispspace;

			uniform float seed = 0;

			normal Nf  = normalize(ntransform(sys, N));
		    point  PP  = noise_freq * transform (sys, P);	
			point PPS = spark_freq * transform(sys, P);
			point Pp = transform(sys, P);
			vector dQu = vtransform(sys, dPdu*du*disp_freq);
			vector dQv = vtransform(sys, dPdv*dv*disp_freq);
			/* sparks */
			float fw = fwidth_mult * filterwidthp(PPS);
			float spark_noise = spark_noise(PPS, spark_amp, octaves, gain, lacunarity, fw);
			spark_noise = smoothstep(0.2, 0.7, spark_noise);
			float angle_mask = 1 - K_angle_mask * (-In) . Nf;
			angle_mask = 1;
			float noise_mask = noise (PPS, seed);
			noise_mask = smoothstep(0.2, 0.8, noise_mask);
			float sparks = angle_mask * noise_mask * spark_noise;
			sparks *= pockmarks; // scale it down
			sparks *= -1; // invert for pockmarks
			/* end sparks */
			float  filterwidth = filterwidthp(PP);
			float  turbf = turbulence(PP, filterwidth, octaves, lacunarity, gain);
			turbf *= snow_noise_scale;

			pxr_brownian(disp_freq, octaves,lacunarity, gain, 
						dispscale, 1, Pp, dQu, dQv, hump);

			getSnow2(Nn, snow_ht, turbf, snowmask);

			turbf += sparks;

			float finalval;

			if(displace_all != 0) { // ignore drift?
				finalval = Km*turbf;
			} else {
				finalval = mix(non_snow_disp*sparks, Km*turbf, snowmask);
			}

			P += normalize(N) * finalval;
			N = normalize(calculatenormal(P));

		}

		public void surface(output color Ci, Oi) {
			float scos, ssin, dot_i, maxrim, nhits =0;
			vector lightdir;
			color texcol = 1, colcomp = 1, dirtcol = 1,
			specweight = 1, shiny = 0, hitc = 0, rt = color(0);
            color inshadow, noshadow;
		    vector In = normalize(I);	
			vector V = -In;
			normal Nf = faceforward(normalize(N), I);
			normal Ns = shadingnormal(N);
            vector raydir = 0;
			vector Ln, H, D;
			if(dispspace == "")
				sys = "shader";
			else
				sys = dispspace;
		    varying point  PP  = spark_freq * transform (sys, P);	
            // we send R to environment light
			vector R = reflect(I,normalize(N));
            // construct ioint data for ptex
            point facedata;
            facedata[0] = u;
            facedata[1] = v;
            facedata[2] = __faceindex;
            normal imageN = normal(0,0,1);
            float Ei;
            float nondiffuse, nonspecular;
            color specK=0, C=0;
            float opacity = 1;

            float ndotv = Ns.V;
            float NH, NHSQ, NH2, Dd, Gg, VH, LN, Ff, tmp;

            ///////////////////////////////////////////////////////////////////////////
            /* shadow pass */
            ///////////////////////////////////////////////////////////////////////////
            illuminance("shad", P) {
            		varying color _shadow;
                    varying float _atten;
                    varying color o_atten = 1;
            		uniform float li_disableshadpass;
            	    // import shadow	
            		lightsource("_shadow", _shadow);
                    // bring in attenuation without any decay
                    lightsource("_atten", _atten);
            		lightsource("disableshadpass", li_disableshadpass);
    
            		if(li_disableshadpass != 1){
            			o_shadow = (_shadow[0]+_shadow[1]+_shadow[2])/3;
                        o_shadow = mix(color(0), color(1), o_shadow);
                        //take attenuation and clamp it
                        o_atten = clamp(o_atten, color(0), color(1));
                        //incorporate attenuation with shadow
                        o_shadow *= o_atten;
            		} else {
            				o_shadow = 1;
            		}
            }

            ///////////////////////////////////////////////////////////////////////////
			/* texture map */
            ///////////////////////////////////////////////////////////////////////////
            float ss = repeat(s, srep);
            float tt = repeat(t, trep);
            filterregion fr;
            fr->calculate2d(ss,tt);
            fr->scale(blurScale);

            if(mapname != "") {
               opacity = texture(mapname[3], s, t);
			   texcol = texture(mapname, fr);
               texcol = mix(texcol, color(1), tmapw);
               //Oi = Os * (opacity - 1);
               if(alphacomp != 0) {
               texcol = mix(texcol, background, 1-opacity);
               }
			}

			if(dirtmap != "") {
				dirtcol = texture(dirtmap, fr);
				dirtcol = 1 - dirtcol; // invert it
			}

            //if(maketexlin != 0) { texcol = sRGB_decode(texcol); }
            if(maketexlin != 0) { texcol = cpow(texcol, 1/0.454545); }

            Oi = Os;
            if(imagetransparency != 0) {
               Oi = Os * mix(color(1), color(0), 1-opacity);
            }
			// SPARKS
			uniform float seed = 0;

			float fw = fwidth_mult * filterwidthp(PP);

			float spark_noise = spark_noise(PP, spark_amp, spark_octaves, spark_gain, spark_lacunarity, fw);
			spark_noise = smoothstep(0.05, 0.5, spark_noise);

			float angle_mask = 1 - K_angle_mask * (-In) . Nf;
			angle_mask = 1;

			float noise_mask = noise (PP, seed);
			noise_mask = smoothstep(0.1, 0.7, noise_mask);

			float sparks = spark_pow * (angle_mask * noise_mask * spark_noise);
			// END OF SPARKS

            ///////////////////////////////////////////////////////////////////////////
            /* bump */
            ///////////////////////////////////////////////////////////////////////////
            float lum = 0;
            color bumps = 1;
            Nn = normalize(N);
            normal Norig = faceforward(normalize(Ng), I);
            normal Nbump = 0;
            normal Nspec = Nn;

            if(texbump != "")
            {
               bumps = texture(texbump, fr);
               if(makebumplin != 0)
               {
                  bumps = sRGB_decode(bumps);
               }
                lum = bumps[0]*0.3 + bumps[1]*0.6 + bumps[2]*0.1;
             	Nn = calculatenormal(P + normalize(N) * lum * Kb);
             	Nbump = faceforward(normalize(Nn), I);
             	Nspec = normalize(mix(Norig,Nbump, clamp(Kb, 0.0, 1.0)));
            }

			if(bumpnoise != 0) {
             	Nn = calculatenormal(P + normalize(N) * clamp(sparks * Kb,0.0,1.0));
             	Nbump = faceforward(normalize(Nn), I);
             	Nspec = normalize(mix(Norig,Nbump, clamp(Kb, 0.0, 1.0)));
			}

			color sparks_out = spark_spec * sparks * category_specular(Nf, -In, 
					spark_roughness);
			filterregion kfr;
			kfr->calculate3d(transform("shader", P));
			float kernel = knoise(kfr, 64, "octaves", 4);
			color sparkle_filter = mix(sparkle1, sparkle2, kernel);
			sparks_out *= sparkle_filter;
			//o_spark = sparks_out;
			o_spark = sparks_out;

			color add_spec = Ks_add * add_spec_color * category_specular(Ns, -In, 
					add_spec_roughness);
			o_add_spec = clamp(add_spec, color(0), color(1));

            ///////////////////////////////////////////////////////////////////////////
			/* diffuse illumination */
            ///////////////////////////////////////////////////////////////////////////
			color k = color(0);
			color di = color(0);
			color cldiff = color(0);
			o_diffusedirect = color(0);

			illuminance("-environment", P, Ns, PI/2) {
				float nondiffuse;
				if(0 == lightsource("__nondiffuse", nondiffuse))
					nondiffuse = 0;
				if(nondiffuse < 1)
				{
					k = Kd * (1-nondiffuse) * Nspec.normalize(L);
					if(0 == lightsource("_shadow", inshadow))
						inshadow = 0;
					if(0 == lightsource("_noshadow", noshadow))
						noshadow = Cl;
					cldiff = noshadow - Cl;
					o_diffusedirect += k*noshadow;
					di += k*cldiff;
				}
			}

            /////////////////////////////////////////////////////////////////////////////
            // indirect illum
            ///////////////////////////////////////////////////////////////////////////
            color indirect = color(0);
            color refl_occ = color(0);
            color colorbleed = color(0);
            color sss = color(0);
            color occl = color(0);
            color kd = color(0);
            o_colorbleeding = 0;
            o_diffusecolor = 0;
            o_refl_occ = 0;
            o_subsurface = 0;
            o_occlusion = 1;
            varying color custom_albedo = clamp(texcol, color(0), color(0.999));
            illuminance("environment", P,
                        "lightcache", "refresh",
                        "send:light:nvec", Ns,
                        "send:light:refl", R,
                        "send:light:facedata", facedata,
                        "send:light:imageN", imageN,
                        "send:light:myalbedo", custom_albedo)
            {
               	lightsource("o_refl_occ", refl_occ);
                  o_refl_occ = clamp(refl_occ, color(0), color(1));
               	lightsource("o_colorbleeding", colorbleed);
                  o_colorbleeding = clamp(colorbleed, color(0), color(1));
               	lightsource("o_subsurface", sss);
                  o_subsurface = clamp(sss, color(0), color(1));
               	lightsource("o_occlusion", occl);
                  //o_occlusion = clamp(occl, color(0), color(1));
               	indirect += Cl;

				if(nonenv == 1) {
					o_refl_occ = 0;
					o_colorbleeding = 0;
					o_subsurface = 0;
					indirect = 0;
				}
            }
            ///////////////////////////////////////////////////////////////////////////
            /* sum up illumination */
            ///////////////////////////////////////////////////////////////////////////

			float non_snowmask = 1 - snowmask; // get black in the snow areas
            kd = Kd * (Cs * texcol * diffusecolor); // sum all sources of color
            o_colorbleeding *= kd;
            o_subsurface *= kd;
            o_refl_occ *= kd;
            o_occlusion *= kd;
            color col = (o_diffusedirect - di) + indirect;
			color illumsurfcolor;
			illumsurfcolor = (Ka * ambient()) + (kd*col)*dirtcol;
            o_diffusecolor = kd;

            ///////////////////////////////////////////////////////////////////////////
            /*Specular */
            ///////////////////////////////////////////////////////////////////////////

        	// Compute reflection roughness values.
			float kss;
        	float refl_roughness_u = roughness(specular_gloss);
        	float refl_roughness_v = refl_roughness_u * max(anisotropy, 0.01);
        
        	if (refl_roughness_u >= 80.0)
        		refl_roughness_u = 80.0 + sqrt(refl_roughness_u - 80.0);
        	if (refl_roughness_v >= 80.0)
        		refl_roughness_v = 80.0 + sqrt(refl_roughness_v - 80.0);

            color highlights_component = 0;
            color hcolor = (is_metallic == 1) ?
                kd*hilightcolor : hilightcolor; 

           if(reflectionSpecularity > 0 && specmode == 0) {
               float computedEcc = eccentricity * eccentricity - 1;
               if(computedEcc > -.00001 && computedEcc < .00001)
                   computedEcc = .00001;
               computedEcc = 1.0 / computedEcc;
               illuminance("-environment", P, Ns, PI/2)
               {
                  if(0 == lightsource("__nonspecular", nonspecular))
                     nonspecular = 0;
                  if(nonspecular < 1)
                  {
                     color noshadow;
                     if(0 == lightsource("_cl_noshadow", noshadow))
                        noshadow = Cl;
                     rfmBlinnSpecular(L, Cl, noshadow, 1-nonspecular,
                             Nspec, V, ndotv, computedEcc, rolloff,
                             reflectionSpecularity, m_shadingCtx->m_raydepth,C, specK);
                     C *= hcolor;
                  }
              }
           } else if(reflectionSpecularity >0 && specmode == 1) {
               C = hcolor * specular_highlight(In,Nspec,V, refl_roughness_u, refl_roughness_v);
               specK = 0;
           } else if(reflectionSpecularity >0 && specmode == 2) { // gloss2
			   illuminance("-environment", P, Ns, PI/2)
			   {
			   		vector H = normalize(normalize(L)+V);
			   		kss = ( smoothstep(.72-w, .72+w, pow(max(0,Nspec.H), 1/specular_gloss)));
			   		C += hcolor * (Cl * kss);
			   		specK = 0;
			   }
		   }
                        
            ///////////////////////////////////////////////////////////////////////////
			/* spec map */
            ///////////////////////////////////////////////////////////////////////////
			if(specmap != "") {
			  specweight = texture(specmap, fr);
              specweight = mix(specweight, color(1), smapw);

              if(makespeclin != 0)
              {
                    specweight = sRGB_decode(specweight);
              }
			} 

            //C = clamp(C, color(0), color(1));
			o_specular = (Ks * specweight) * C;
            o_specularshadow = (Ks * specweight) * specK;

            ///////////////////////////////////////////////////////////////////////////
            /* reflection */
            ///////////////////////////////////////////////////////////////////////////
            vector Rn, T0, T1;
            vector En, T2, T3;
	        color Cr = 0, rho = 0;
            normal myN = normalize(N);

			// grazing angles for snow reflections
			float grangle = Nn.V;
			grangle = clamp(1-grangle, 0, 1);
			grangle = pow(grangle, graze);
			grangle *= frnl; // overall weight

            EnvReflFrame(In, myN, angle, yup_to_zup, envmapspace, En, T2, T3);
            ReflFrame(In, myN, angle, yup_to_zup, envmapspace, Rn, T0, T1);

			// fresnel for non-snow parts
            rho = FrnlBlend(kd, frnl_blend, abs(In.myN));
			rho *= frnl; // overall weight

            if(Kr > 0 && env_ctx == "all" ) {
                rt += GatherRefl(P, Rn, T0, T1, envmap, yup_to_zup, envmapspace, samples, fadereflections, fadescale, fadeexp, 0 );
            } 
             else if(Kr > 0 && env_ctx == "env")  {
	             rt = EnvRefl(En, T2, T3, envmap);
            } 
            else if(Kr > 0 && env_ctx == "surface") {
                rt += GatherRefl(P, Rn, T0, T1, envmap, yup_to_zup, envmapspace, samples, fadereflections, fadescale, fadeexp, 1 );
            } else {
                rt = 0;
            }

            color reflhsv = 0;
            color reflrgb = 0;
            color myrt = rt;
            reflhsv = ctransform("hsv", myrt);

            if(refltint != 0) { // desat reflections
                setcomp(reflhsv, 1, refl_saturation);
                reflrgb = ctransform("hsv", "rgb", reflhsv);
                rt = reflrgb;
                rt = mix(rt, kd, refl_blend);
            }

            /* reflection map */ 
            color Rmap = 1;
            if(reflmap != "")
            {
               Rmap = texture(reflmap, fr);
               Rmap = mix(Rmap, color(1), rmapw);
               // linear map
               if(makerefllin != 0)
               {
                  Rmap = sRGB_decode(Rmap);
               }
            }
			// we need the snow to have different reflection values
			// invert the mask so it's black in snow areas, and thus negates the refl map
			float rmap_fl;
			rmap_fl = max(Rmap[0],snowmask); // the reflection map has no effect on snow portions
            o_reflection = Kr * (rmap_fl * rt);

				
            ///////////////////////////////////////////////////////////////////////////
			/* Rim light */
            ///////////////////////////////////////////////////////////////////////////
			illuminance("-environment", P, Nf, PI) {
			    scos = Nf.V;
			    ssin = sqrt(1-pow(scos,2));
                color noshadow;
                if(0 == lightsource("_cl_noshadow", noshadow))
                noshadow = Cl;
                color cldiff = noshadow - Cl;
				Ci = noshadow * lightcolor;
				lightsource("__nonspecular", nonspecular);
				if(nonspecular < 1) {
				Ci *= (1 - nonspecular);
				
				Ln = normalize(L);
				scos = Nf.V;
				ssin = sqrt(1-pow(scos,2));
									
				shiny += pow(ssin, (1/rimstrength*5))*(Ln.Nf) * noshadow * rimcolor;					
				shiny += pow(scos, (rimstrength*5)) * (Ln.Nf) * noshadow * color(0);
				}
			}
			o_rimcolor = Krim*shiny;


			// snow mask
			// getSnow(Nn, snow_ht, snow_exp, driftDir, snowmask);
			// snow gets Buratti shading model
			color Buratti_C = snowcolor * Buratti(Nn, In, snowK, sden) + indirect;

			illumsurfcolor =  mix(illumsurfcolor, (Buratti_C + sparks_out + add_spec + o_rimcolor), snowmask);
			o_reflection   =  mix(rho*o_reflection, grangle*o_reflection, snowmask);
			//o_reflection *= rho;


			/* SSS */
			color sssval = 0;
			color s_albedo = 0;
			s_albedo = mix(color(0), albedo, snowmask);
			if(sssmap != "") {
				sssval = subsurface(P, Ns, "filename", sssmap,
								"albedo", s_albedo, "diffusemeanfreepath", dmfp,
								"ior", ior, "unitlength", unitlength,
								"smooth", smooth,
								"maxsolidangle", maxsolidangle);
				o_subsurface = sssval;
				sssval *= snowmask;
				illumsurfcolor += sssval;
				//finalsum = clamp(sssval + finalsum, color(0), color(1));
				//illumsurfcolor = mix(illumsurfcolor, (sssval + sparks_out + add_spec), snowmask);
			}

            /* Final Values */
            color finalsum = ((illumsurfcolor + o_reflection + o_rimcolor) + o_specular);

            ///////////////////////////////////////////////////////////////////////////
            /* bakeing & reading pointclouds  */
            ///////////////////////////////////////////////////////////////////////////
            color bakecolor = 0;
            normal bk = (bakenormal == 0) ? normal(0,0,0) : Nspec;
            point imageP = (u,v,0);
            // normal output
        	float nmapx = 0.5*(1+xcomp(Nspec));
        	float nmapy = 0.5*(1+ycomp(Nspec));
        	float nmapz = 0.5*(1+zcomp(Nspec));
			float a = area(P, "dicing");
        	o_normals = color(nmapx, nmapy, nmapz);

            if(imagebake != 0) {
                imageP = point(s,t,0);
            }
            if(m_shadingCtx->m_raydepth < 1 && bakemap != "") {
                if (bakedata == "occlusion")
                {
                    bakecolor = o_occlusion;
                } else if (bakedata == "illumination")
                {
                    bakecolor = finalsum;
                } else if (bakedata == "colorbleeding")
                {
                    bakecolor = o_colorbleeding;
                } else if (bakedata == "refl_occ")
                {
                    bakecolor = o_refl_occ;
                } else if (bakedata == "normals")
                {
                    bakecolor = o_normals;
                }

                if(bakemode == "ptex") 
				{
                   bake3d(bakemap, displaychannel, facedata, imageN, "coordsystem",
                        "_disable", "interpolate", 1,
                        displaychannel, bakecolor);

                } 
				else if (bakemode == "uv")
				{
	               bake3d(bakemap, displaychannel, imageP, normal(0), "interpolate", 
	               		interpolate, displaychannel, bakecolor);

                }
				else if (bakemode == "ptc" && bakedata == "_area,_radiosity,Cs")
				{
                    bake3d(bakemap, displaychannel, P, Nn, "interpolate", 1,
                          "_area", a, "_radiosity", Ci, "Cs", finalsum);
                }
				else if (bakemode == "ptc")
				{
					bake3d(bakemap, displaychannel, P, bk, "interpolate", interpolate,
							displaychannel, bakecolor);
				}
            }

			// spherical harmonics
			uniform color envmapshcoeffs[];
			float dirvisshcoeffs[];
			float incradshcoeffs[];
			color envcol = 0;
			color irradcol = 0;
			float ok;

            filterregion ptxfr;
            ptxfr->calculate2d(u,v,du,0,0,dv);
            color ptxdata = 0;
            if(ptexmap != "") {
               ptxdata = ptexture(ptexmap, 0, __faceindex, ptxfr); 
			   finalsum = ptxdata;
               //Ci = Oi * ptxdata;
            } else if (ptcmap != "") { 
			   textureinfo(envmap, "shcoeffs", envmapshcoeffs);
               ok = texture3d(ptcmap, P, Nn, "_dirvisshcoeffs", dirvisshcoeffs);
			   if(ok != 0) {
				   stdrsl_SphericalHarmonic envmapsh, dirvissh, radvissh;
				   envmapsh->createFromArray(envmapshcoeffs);
				   dirvissh->createFromArray(dirvisshcoeffs);
				   radvissh->createFromArray(incradshcoeffs);
				   envcol = envmapsh->convolve(dirvissh, 0);
				   irradcol = envmapsh->convolve(radvissh, 0);
				   o_occlusion = envcol;
				   o_colorbleeding = irradcol;
			   }
            } 
            /* final output*/
		    Ci = finalsum * Os;
			Oi = Os;
		    //Ci = Oi * ((illumsurfcolor + o_reflection) + (o_specular + o_rimcolor));

		}
}

	
