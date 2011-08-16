//
// aw_prmanEnvLight.sl - Environment light for indirectdiffuse,
// occlusion, reflection occlusion, caustics and sss
// Alan Warren, 2011.
//

#include <stdrsl/ShadingContext.h> // for ShadingUtils
#include "aw/aw_PrmanEnvLight.h"
#include "normals.h"


class aw_prmanEnvLight( uniform string filename = "";
                        uniform float doall = 0;
                        varying float maxdist = 1e15; 
		                varying float falloff = 0;
		                varying float falloffmode = 0;
		                varying float samplebase = 0;
		                varying float bias = 0.009;
		                varying float clamp = 1;
		                varying float maxsolidangle = 0.05;
		                varying float maxvariation = 0;
		                uniform string envmap = "";
		                uniform string hitsides = "both";
		                varying float outswitch = 1;
		                uniform string distribution = "cosine";
		                uniform float sortbleeding = 1;
		                varying float envswitch = 0;
		                varying color albedo = color(0.8, 0.79, 0.75);
		                varying color dmfp = color(8.51, 5.57, 3.95);
		                varying float ior = 1.5;
		                varying float unitlength = 1;
		                varying float smooth = 0;
		                uniform float samples = 0;
		                varying float causticSW = 1;
		                uniform string causticfile = "";
		                uniform float estimator = 100;
		                uniform float coneangle = 0.200000003;
		                uniform float albedo_switch = 0;
		                uniform float rasterresolution = 12;
		                uniform string myfilename = "";
		                varying float bake_file_index = 0;
                        uniform float glossy_refl = 0;
		                uniform string envspace = "world";
		                varying color filter = color(1, 1, 1);
		                uniform string __category = "environment";
                        )

{
    varying normal m_Ns;
    constant float maxddepth; 
    constant float sides;
    string the_ptc;
    uniform float raydepth;
    uniform float ddepth;
    varying color envcol;

    public void construct()
    {
        attribute("trace:maxdiffusedepth", maxddepth);
        attribute("Sides", sides);
    }

    public void begin()
    {
        if(bake_file_index == 0)
            the_ptc = filename;
        else
            the_ptc = myfilename;

        rayinfo("depth", raydepth);
        rayinfo("diffusedepth", ddepth);

        normal nn, nf;
        vector in;
        stdrsl_ShadingUtils sutils;
        sutils->GetShadingVars(raydepth, sides, nn, in, nf, m_Ns);
    }

    public color aw_occlusion()
    {
            // environment occlusion
            color occl = 0;
            vector illumdir = (0,0,0);
            color envresult = 0;
            vector envdir = 0;
            if(envswitch == 0 && outswitch == 0 || doall == 1 && envswitch == 0) {
                occl = occlusion(Ps, m_Ns, 0, "pointbased", 1,
                        "filename", the_ptc, "hitsides", hitsides,
                        "distribution", distribution,
                        "maxdist", maxdist, "falloff", falloff,
                        "falloffmode", falloffmode,
                        "samplebase", samplebase, "bias", bias,
                        "clamp", clamp, "maxsolidangle", maxsolidangle,
                        "maxvariation", maxvariation,
                        "rasterresolution", rasterresolution,
                        "illuminationdir", illumdir,
                        "environmentmap", envmap,
                        "environmentspace", envspace,
                        "environmentcolor", envcol,
                        "environmentdir", envdir);
                // sum
                envresult = envcol * (1 - occl); 
                return filter * envresult;
                //lr->occl = product;
                //lr->illuminationdir = illumdir;

            // standard occlusion
            } else if(envswitch == 1 && outswitch == 0 || doall == 1 && envswitch == 1) {
                occl = occlusion(Ps, m_Ns, 0, "pointbased", 1,
                        "filename", the_ptc, "hitsides", hitsides,
                        "distribution", distribution,
                        "maxdist", maxdist, "falloff", falloff,
                        "falloffmode", falloffmode,
                        "samplebase", samplebase, "bias", bias,
                        "clamp", clamp, "maxsolidangle", maxsolidangle,
                        "maxvariation", maxvariation,
                        "rasterresolution", rasterresolution);
                // sum
                return filter * (1 - occl);

            }
    }

    public color aw_indirectdiff()
    {
        // indirect diffuse env
        color indirectdiff;
        vector illumdir = (0,0,0);
        if(envswitch == 0 && outswitch == 1 || doall == 1 && envswitch == 0) {
            indirectdiff = indirectdiffuse(Ps, m_Ns, 0, 
                    "pointbased", 1,
                    "filename", the_ptc,
                    "hitsides", hitsides,
                    "colorhitsides", "front",
                    "clamp", clamp, "sortbleeding", sortbleeding,
                    "maxdist", maxdist, "falloff", falloff,
                    "falloffmode", falloffmode,
                    "samplebase", samplebase, "bias", bias,
                    "maxsolidangle", maxsolidangle,
                    "rasterresolution", rasterresolution,
                    "illuminationdir", illumdir,
                    "maxvariation", maxvariation,
                    "environmentcolor", envcol,
                    "environmentspace", envspace,
                    "environmentmap", envmap);
            // sum
            //id = lr->illuminationdir;
            return filter * indirectdiff;
            //return filter * envcol;
            //lr->indirectdiff = product;
            //lr->illuminationdir = illumdir;

        // indirect diffuse no env
        } else if(envswitch == 1 && outswitch == 1 || doall == 1 && envswitch == 1) {
            indirectdiff = indirectdiffuse(Ps, m_Ns, 0, 
                    "pointbased", 1,
                    "filename", the_ptc,
                    "hitsides", hitsides,
                    "colorhitsides", "front",
                    "clamp", clamp, "sortbleeding", sortbleeding,
                    "maxdist", maxdist, "falloff", falloff,
                    "falloffmode", falloffmode,
                    "samplebase", samplebase, "bias", bias,
                    "rasterresolution", rasterresolution,
                    "maxsolidangle", maxsolidangle,
                    "maxvariation", maxvariation);
            // sum
            return filter * indirectdiff;
            //lr->indirectdiff = product;
        }
    }

    public color aw_sss(point pp)
    {
        color albedoresult = 0;
        if(albedo_switch == 0)
            albedoresult = albedo;
        else
            //albedoresult = myalbedo;
            //TODO: use co-shader for user supplied value
            albedoresult = albedo;

        // sss
        color subsurf = 0;
        if(outswitch == 2) {
            subsurf = subsurface(pp, m_Ns, "filename", the_ptc,
                    "albedo", albedoresult, "diffusemeanfreepath", dmfp,
                    "ior", ior, "unitlength", unitlength,
                    "smooth", smooth, "maxsolidangle", maxsolidangle);
            // sum
            return filter * subsurf;
            //lr->subsurf = product;
        } 
            
    }

    public color aw_photons()
    {
        color irrad = 0;
        color radiosity = 0;
        if(outswitch == 3 && causticSW == 0) {
            irrad = filter * photonmap(causticfile, Ps, m_Ns, 
                    "estimator",estimator);
             return irrad;
        // irradiance photons
        } else if(outswitch == 3 && causticSW == 1) {
            if(ddepth == maxddepth) { // brick map
                texture3d(the_ptc, Ps, m_Ns, "_radiosity", radiosity);
                return radiosity;
            } else {
                return filter * indirectdiffuse(Ps, m_Ns, samples, "samplebase", 0.1,
                    "maxvariation", maxvariation);
            }
        }
    }

    public color aw_refl_occ(vector dir) 
    {
        // reflection occlusion
        color occlude = 0;
        color envresult = 0;
        vector illumdir = (0,0,0);
        vector envdir = 0;
        if(outswitch == 4 && envswitch == 0 || doall == 1 && envswitch == 0) {
            occlude = occlusion(Ps, dir, 0, "pointbased", 1,
                    "filename", the_ptc, "hitsides", hitsides,
                    "distribution", distribution,
                    "maxdist", maxdist, "falloff", falloff,
                    "falloffmode", falloffmode,
                    "samplebase", samplebase, "bias", bias,
                    "clamp", clamp, "maxsolidangle", maxsolidangle,
                    "maxvariation", maxvariation,
                    "rasterresolution", rasterresolution,
                    "illuminationdir", illumdir,
                    "environmentmap", envmap,
                    "environmentspace", envspace,
                    "environmentcolor", envcol,
                    "environmentdir", envdir);
            // sum
            envresult = envcol * (1 - occlude); 
            return filter * envresult;
            //lr->refl_occ = product;
            //lr->illuminationdir = illumdir;

        } else if(outswitch == 4 && envswitch == 1 || doall == 1 && envswitch == 1) {
            occlude = occlusion(Ps, dir, 0, "pointbased", 1,
                    "filename", the_ptc, "hitsides", hitsides,
                    "distribution", distribution,
                    "maxdist", maxdist, "falloff", falloff,
                    "falloffmode", falloffmode,
                    "samplebase", samplebase, "bias", bias,
                    "clamp", clamp, "maxsolidangle", maxsolidangle,
                    "maxvariation", maxvariation,
                    "rasterresolution", rasterresolution);
            // sum return filter * (1 - occlude);

        }
    }

    public color aw_refl_indirect(vector dir)
    {
        // glossy refl env
        color indirectdiff;
        color envresult = 0;
        vector illumdir = (0,0,0);
        if(envswitch == 0 && outswitch == 1 && glossy_refl == 1) {
            indirectdiff = indirectdiffuse(Ps, dir, 0, 
                    "pointbased", 1,
                    "filename", the_ptc,
                    "hitsides", hitsides,
                    "distribution", distribution,
                    "coneangle", radians(coneangle),
                    "colorhitsides", "front",
                    "clamp", clamp, "sortbleeding", sortbleeding,
                    "maxdist", maxdist, "falloff", falloff,
                    "falloffmode", falloffmode,
                    "samplebase", samplebase, "bias", bias,
                    "maxsolidangle", maxsolidangle,
                    "rasterresolution", rasterresolution,
                    "illuminationdir", illumdir,
                    "maxvariation", maxvariation,
                    "environmentspace", envspace,
                    "environmentcolor", envcol,
                    "environmentmap", envmap);
            //envresult = envcol * (1 - indirectdiff);
            // sum
            return filter * (indirectdiff);
            //return filter * envresult;

        // glossy refl no env
        } else if(envswitch == 1 && outswitch == 1 && glossy_refl == 1) {
            indirectdiff = indirectdiffuse(Ps, dir, 0, 
                    "pointbased", 1,
                    "filename", the_ptc,
                    "hitsides", hitsides,
                    "distribution", distribution,
                    "coneangle", radians(coneangle),
                    "colorhitsides", "front",
                    "clamp", clamp, "sortbleeding", sortbleeding,
                    "maxdist", maxdist, "falloff", falloff,
                    "falloffmode", falloffmode,
                    "samplebase", samplebase, "bias", bias,
                    "rasterresolution", rasterresolution,
                    "maxsolidangle", maxsolidangle,
                    "maxvariation", maxvariation);
            // sum
            return filter * indirectdiff;
            //lr->indirectdiff = product;
        }
    }


    public void light( output vector L; output color Cl)
    {
       L = (0,0,0);
       Cl = (0,0,0);
    }

}// aw_prmanEnvLight
