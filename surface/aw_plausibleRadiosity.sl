/************************************************************************
 * aw_plausibleRadiosity.sl - shader for baking radiosity. This shader
 * is injected in the RIB stream via an Ri Filter by my envlight HDA
 *
 * Author: Alan Warren (bluemoonshine@gmail.com)
 *
 * $Revision: 1.1 $    $Date: 2008/10/20 $
 *
 ************************************************************************/


#include <stdrsl/ShadingContext.h>
#include "aw/aw_TexturingUtils.h"


surface
   aw_plausibleRadiosity(string bakefile = "";
                    string displaychannels = "_area,_radiosity,Cs,Os";
                    string texturename = ""; 
                    string bakeshader = "radiositytexture";
                    float Ka = 1, Kd = 1;
                    float __faceindex = 0;
                    float houdini = 1;
                    float texblur = 1;
                    output uniform string mytexture = "")
{
            stdrsl_ShadingContext m_shadingCtx;
            aw_TexturingUtils a_tex;
            a_tex->aw_initFr(texblur, houdini, __faceindex);
            m_shadingCtx->init();
            color irrad, tex = 1, diffcolor;
            normal Nn = normalize(N);
            float a = area(P, "dicing");
            color c = 1;
            string texmap = texturename;

                        if(bakeshader != "")
                        {
                                shader brt = getshader(bakeshader);
                                if(brt != null) {
                                    brt->getTexture(a_tex, c);
                                }
                        } 
                        else if(texmap != "")
                        {
                            c = color texture(texturename, a_tex->m_fr);
                            //tex = texture(texmap, s, t, "filter", "gaussian",
                            //"width", 0.75, "lerp", 1);
                        }
            
            irrad = Ka*ambient() + Kd*diffuse(Nn);
      
            diffcolor = Cs * c;
      
            Ci = irrad * Os;
            Oi = Os;
      
            if (a > 0) 
                        {
              bake3d(bakefile, displaychannels, P, Nn, "interpolate", 1,
                     "_area", a, "_radiosity", Ci, "Cs", diffcolor, "Os", Oi);
                        }
} 



