#ifndef aw_TexturingUtils_h
#define aw_TexturingUtils_h

#include <stdrsl/Math.h>

struct aw_ShadingUtils
{
    public void aw_GetShadingVars(float textureBlur;
                                  float houdini;
                                  output float adu;
                                  output float adv;
                                  output float au; 
                                  output float av;
                                  output float as;
                                  output float at;
                                  output filterregion fr;
                                  output filterregion ptxfr;
                                  output point facedata;
                                  output float __faceindex;)
    {
        extern float du;
        extern float dv;
        extern float u;
        extern float v;
        extern float s;
        extern float t;
        adu = du;
        adv = dv;
        au = u;
        av = v;
        as = s;
        at = t;
        facedata[0] = u;
        facedata[1] = v;
        facedata[2] = __faceindex;
        ptxfr->calculate2d(au,av,adu,0,0,adv);
        fr->calculate2d(as, (houdini>0)?(1-at):at); 
        fr->scale(textureBlur);
    }
};
//
// aw_ShadingContext:
//   a container for geometric and other context relevent to BSDFs
//
struct aw_TexturingUtils
{
    varying float m_adu = 0;
    varying float m_adv = 0;
    varying float m_au = 0;
    varying float m_av = 0;
    varying float m_as = 0;
    varying float m_at = 0;
    varying filterregion m_fr;
    varying filterregion m_ptxfr;
    varying point m_facedata = 0;
    varying float m__faceindex = 0;

    public void aw_initFr( float blurScale;
                           float houdini;
                           float __faceindex;)
                           //filterregion fr;
                           //filterregion ptxfr;)
    {
        aw_ShadingUtils asutil;
        asutil->aw_GetShadingVars(blurScale, houdini, m_adu, m_adv, m_au, m_av, 
                m_as, m_at, m_fr, m_ptxfr, m_facedata, m__faceindex);

    }


    public void aw_displace(vector input; float mag; color tex; string mode;) 
    {
        extern point P;
        extern normal Ng;
        extern normal N;

        point tmpP = P; 
        normal Nn  = normalize(N);
        normal Ngn = normalize(Ng);
        normal normalDelta = Nn - Ngn; //Normalized difference for polys
        float referenceMag = length(vector "shader" (1,0,0));
        float level = tex[0];
        float strength = ((level * 2) - 1) * mag;
        P += input * strength * referenceMag; 
        N = normalize(calculatenormal(P)) + normalDelta;
        if(mode == "bump") P = tmpP;
    }

    public void aw_ptexDisplace(vector input; float mag; color ptx; string mode;)
    {
        extern point P;
        extern normal Ng;
        extern normal N;

        point tmpP = P;
        point objP = transform("object", P);
        normal Nn  = normalize(N);
        normal Ngn = normalize(Ng);
        normal normalDelta = Nn - Ngn;
        float referenceMag = length(vector "shader" (1,0,0));
        objP = objP + (mag * normal(ptx[0],ptx[1],ptx[2]) * referenceMag);
        P = transform("object", "current", objP);
        N = normalize(calculatenormal(P)) + normalDelta;
        if(mode == "bump") P = tmpP;
    }

};

#endif
