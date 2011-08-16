/************************************************************************
 * aw_glass_v5.sl - Glass shader with absorption, glossy spec and color shadows
 *
 * Author: Alan Warren (bluemoonshine@gmail.com)
 *
 *	Credit: Attenuation ported from Mario Marengo's VEX Glass on od[force] except the way 
 *  in which ray length is obtained, which was an idea posted by Pixar's staff in the 
 *  prman forums. Also, Pixar's glassrefr.sl was used for ideas. In other words I didn't 
 *  do much but I can take credit for hacking with the opacity at the end for 
 *  transmission shadows :) 
 *
 * $Revision: 1.1 $    $Date: 2008/09/14 01:59:46 $
 *
 ************************************************************************/

#include "../include/rmannotes.h"
#include"../include/normals.h"


class aw_glass_v5( float Kr = 1; 
		           float Kt = 1; 
		           float Ks = 1;
                   color surface_id = 1;
                   float ior = 1.5;
                   float frnl = 1;
                   float rough = 0.5;
		           float refr_blur = 0;
                   float refl_blur = 0;
                   float blurScale = 1;
                   float srep = 1;
                   float trep = 1;
                   float samples = 1;
                   string mapname = "";
                   string envmap = "";
		             /*absorbtion*/
                   float absorption = 100;
		           float fadeexp = 1;
		           float fadeamplitude = 1;
		           color tint = color(0.42,0.87,0.84);
		           /* outputs */
		           output varying color o_reflection = 0;
                   output varying color o_id = 0;
		           output varying color o_transmission = 0;
                   output varying color o_texture = 0;)
{

  public void ReflFrame(vector In; vector Nn; float coneSize;
                        output vector Rn, T0, T1)
  {
     Rn = In - 2 * (In.Nn)*Nn;
     T1 = normalize(Rn ^ Nn);
     T0 = normalize(Rn ^ T1);
     T0 *= coneSize;
     T1 *= coneSize * Nn.Rn;
  }

  private color GatherRefl(point P; vector Rn; vector T0; vector T1; string envmap;)
  {
     color ret = 0;
     color Cret = 0;
     vector rdir = 0;
     uniform float raydepth;
     rayinfo("depth", raydepth);
     uniform float numSamps = (raydepth < 1) ? samples : 1;

     filterregion fr;
     fr->calculate3d(Rn, T0, T1);

     gather("illuminance", P, fr, numSamps, "surface:Ci", Cret, "ray:direction", rdir) {
        ret += Cret;
     } else {
         if(envmap != "") {
            filterregion frEnv;
            frEnv->calculate3d(normalize(rdir));
            ret += environment(envmap, frEnv);
         }
     }
     return ret / numSamps;
  }

  private color FrnlBlend(color vNormal; color vGrazing; float cosTheta)
  {  // schlick approx
     float wt = pow(1-cosTheta, 5);
     return mix(vNormal, vGrazing, color(wt));
  }

				
  public void surface(output color Ci, Oi;) {
		  normal Nn = 0;
          normal n = normalize(N);
		  vector In = normalize(I);
		  normal Nf = faceforward(n, In, n);
		  normal Ns = shadingnormal(N);
		  vector V = -In;   /* view direction */
		  vector R=0, T=0;
          float entering = (In.n <= 0) ? 1 : 0;
		  float eta = (entering == 1) ? (1/ior) : ior;
		  //float entering = (V.Nn > 0.0) ? 1 : 0;
		  float kr, kt, depth;
		  color rt = 0, rhitc = 0, thitc = 0;
		  vector dir = 0;
		  float Rnhits = 0, Tnhits = 0, dt = 0, maxcont = 1;
		  float nonspec = 0;
		  color atten = 1, k = 1, tz = 1;
		  	  
		  uniform float raydepth;
          point hitP = 0;
          float tracedlen = 1.0e30;
		  rayinfo("depth", raydepth);
		  float xsamples = (raydepth == 0) ? samples : 1; 

		  color tinthsv = ctransform("rgb","hsv", tint);
		  color Cabs = ctransform("hsv","rgb", color(tinthsv[0],tinthsv[1],1));	
          color kd = 1;
		  

		  fresnel(In, Ns, eta, kr, kt, R, T);
		  kt = (1 - kr);

		  Oi = 1;
		   
		  /* attenuate refracted rays by the distance they travel
		   * i'm taking the length of the vector from the ray origin
         * to the ray hit point for secondary rays after they've 
         * passed through the surface
         */
		  if(raydepth > 0 && entering == 0) {
            gather("", P, I, 0, 1, "surface:P", hitP)
            {
               tracedlen = length(hitP - P);
            }

            float d = (absorption - tracedlen) / absorption;
            d = smoothstep(0,1,d);
            d = pow(d, fadeexp);

				atten[0] = exp((Cabs[0]-max(Cabs[0],Cabs[1],Cabs[2])-
							min(Cabs[0],Cabs[1],Cabs[2]))*2.0*fadeamplitude*d);
				atten[1] = exp((Cabs[1]-max(Cabs[0],Cabs[1],Cabs[2])-
							min(Cabs[0],Cabs[1],Cabs[2]))*2.0*fadeamplitude*d);
				atten[2] = exp((Cabs[2]-max(Cabs[0],Cabs[1],Cabs[2])-
							min(Cabs[0],Cabs[1],Cabs[2]))*2.0*fadeamplitude*d);
        }
		  
            ///////////////////////////////////////////////////////////////////////////
            /* reflection */
            ///////////////////////////////////////////////////////////////////////////
            vector Rn, T0, T1;
            normal myN = normalize(N);
            ReflFrame(In, myN, refl_blur, Rn, T0, T1);
            color rho = FrnlBlend(kd, color(1), abs(In.myN));
            if(Kr * kr > 0 ) {
                k = mix(color(kr), kr*atten, clamp(Kr,0,1));
                maxcont = max(k[0],k[1],k[2]);
                if(maxcont >0.0) {
                 //o_reflection = k * trace(P, R);
                rho *= frnl;
                rt += rho * GatherRefl(P, Rn, T0, T1, envmap);
                o_reflection = k*rt;
                }
            }

          float ss = repeat(s, srep);
          float tt = repeat(t, trep);
          filterregion fr;
          fr->calculate2d(ss,tt);
          fr->scale(blurScale);
          float texfloat = 1.0;
          color texcol = 1.0;

          if(mapname != "") {
		    texcol = texture(mapname, fr);
          }
          texfloat = (texcol[0] + texcol[1] + texcol[2]) / 3;
          o_texture = texcol;

		 /* Refraction */
         texfloat *= Kt;
		 color ttt = 0;
		 if (Kt * kt > 0) {
          tz = mix(color(kt), kt*atten, clamp(texfloat,0,1));
          maxcont = max(k[0],k[1],k[2]);
          if(maxcont>0.0) {
             //o_transmission = tz * trace(P, T);
		    gather("illuminance", P, T, radians(refr_blur),
		    		xsamples,  "volume:Ci", thitc)
		    {
		    		ttt += thitc;
		    		Tnhits += 1;
		    }
		    ttt /= Tnhits;
		    o_transmission = tz * ttt;
            }
				   
		  }
		  /* generate an opacity value for transmission based shadows. It's color by default will be
		   * the complement of our intended color. To fix this we invert the hue
		   */
		  color opacity = mix(color(kt), kt*atten, clamp(Kt,0,1));
		  color darkmix = mix(color(kr), kr*atten, clamp(Kr,0,1));
		  color shad = tint;
		  //convert to hsv
		  opacity = ctransform("hsv", opacity);
		  darkmix = ctransform("hsv", darkmix);
		  //desaturate completely. 
		  setcomp(opacity, 1, 0);
		  setcomp(darkmix, 1, 0);
		  //convert it back to rgb.
		  opacity = ctransform("hsv", "rgb", opacity);
		  darkmix = ctransform("hsv", "rgb", darkmix);
		  //convert shad to hsv, invert it's hue, then convert it back to rgb
		  color component = ctransform("rgb", "hsv", shad);
		  float reversehue = component[0]+0.5;
		  reversehue = (reversehue > 1) ? reversehue-1.0 : reversehue;
		  setcomp(component, 0, reversehue);
		  shad = ctransform("hsv", "rgb", component);

          // surface id
          o_id = surface_id;
		  
		  uniform string raytype;
		  rayinfo("type", raytype);
		  if(raytype == "transmission") { Oi = (opacity * shad)+darkmix;
           Ci = color(0);
        }
		  
		  Ci =  (Ks*texfloat) * specular(n, V, rough) + (o_reflection+o_transmission);

		 
		}
}
