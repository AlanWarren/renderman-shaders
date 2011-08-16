
/* Shader Utilities */
#define ALPHA 7

//decode from sRGB luma to linear light
float sRGB_decode_f(float f)
{
   float lin;
   if(f <= 0.03928)
      lin = f/12.92;
   else
      lin = pow((f+0.055)/1.055, 2.4);
   return lin;
}

color sRGB_decode(color c)
{
   color d;
   d[0] = sRGB_decode_f(c[0]);
   d[1] = sRGB_decode_f(c[1]);
   d[2] = sRGB_decode_f(c[2]);
   return d;
}

color cpow(color c; float p)
{
    color d;
    d[0] = pow(c[0], p);
    d[1] = pow(c[1], p);
    d[2] = pow(c[2], p);
    return d;
}

/* DspTangentSpace for Vector Displacement
 * credit Philippe Leprince
 */
struct DspTangentSpace
{
    varying vector tx = 0;
    varying vector ty = 0;
    varying vector tz = 0;

    void init( varying point p; varying vector n; )
    {
        // build tangent space
        tz = normalize(n);
        tx = normalize( ( Du(p)*Dv(v) - Dv(p)*Du(v) ) / ( Du(u)*Dv(v) - Dv(u)*Du(v) ) );
        ty = normalize( tz ^ tx );
    }

    vector toTangentSpace( varying vector iv; )
    {
        return vector(tx.iv , tz.iv , ty.iv);
    }

    vector fromTangentSpace( varying vector iv; )
    {
        return vector( iv[0]*tx[0] + iv[1]*tz[0] + iv[2]*ty[0],
                       iv[0]*tx[1] + iv[1]*tz[1] + iv[2]*ty[1],
                       iv[0]*tx[2] + iv[1]*tz[2] + iv[2]*ty[2] );
    }
};

//remap
float remap(float x, a1, b1, a2, b2) {
		return (x*(b2-a2) - a1*b2 + b1*a2) / (b1-a1);
}

color remapc(color x; float a1, b1, a2, b2) {
		float red = (x[0]*(b2-a2) - a1*b2 + b1*a2) / (b1-a1);
		float green = (x[1]*(b2-a2) - a1*b2 + b1*a2) / (b1-a1);
		float blue= (x[2]*(b2-a2) - a1*b2 + b1*a2) / (b1-a1);
		return color(red, green, blue);
}
vector remapv(vector x; float a1, b1, a2, b2) {
		float X = (x[0]*(b2-a2) - a1*b2 + b1*a2) / (b1-a1);
		float Y = (x[1]*(b2-a2) - a1*b2 + b1*a2) / (b1-a1);
		float Z = (x[2]*(b2-a2) - a1*b2 + b1*a2) / (b1-a1);
		return vector(X, Y, Z);
}
		
//color2float
float color2float(color col) {
		return (col[0]+col[1]+col[2])/3;
}

//environment mapping tools
color envfunc(string envname)
{
		color Cr = 0;
		if(envname != "")
		{
				vector Nf = faceforward(normalize(N), I);
				vector V = -normalize(I);
				vector Rray = reflect(-V, Nf);
				Rray = vtransform("shader", Rray);
				Cr = color environment(envname, Rray);
		}
		return Cr;
}
// handy for making +Y the north pole of env lookups
vector yup_to_zup(vector refldir)
{
   return (vector rotate((point refldir), PI/2, point(0,0,0), point(1,0,0)));
}

#define ALPHA 7 
float scaledschlick(float theta, ior, extcoeff;)
{
		float out;
		out = (pow(ior-1, 2)+4*ior*pow(1-theta, 5)+pow(extcoeff, 2)) /
				(pow(ior+1, 2)+pow(extcoeff, 2)); // term
		if(extcoeff > 0 || ior > 1.5)
		out = out - 2 * ior * theta * pow(1-theta, ALPHA);
		
		return out;
}
//float specular BRDF
float specBRDF_f(vector L, N, V; float roughness)
{
   vector H = normalize(L+V);
   return pow(max(0, N.H), 1/roughness);
}
// common raytracing lookup
/* struct ShadingContext {
  uniform float m_raydepth = 0;
  string m_passName = "";
  float m_frameNumber = 0;
  public void InitializeOptions() {
    option("user:passName",m_passName);
    option("Frame",m_frameNumber);
  }
  public void InitializeAttributes() {
    rayinfo("depth",m_raydepth);
  }
 //} */

vector DerivEx(point p; varying float ss, tt;)
{
// consider the 2x2 Jacobian matrix
// J = | Du(ss) Dv(ss) |
// | Du(tt) Dv(tt) |
// compute the determinant of matrix J
float detJ = Du(ss)*Dv(tt) - Du(tt)*Dv(ss);

// inverse(J) = (1/determinant(J)) * | Dv(tt) -Dv(ss) |
// | -Du(tt) Du(ss) |
// [dPdss dPdtt] = [Du(P) Dv(P)] * inverse(J)
//
// this can be translated as:
// dPdss = (Du(P)*Dv(tt) - Dv(P)*Du(tt))/det(J)
// dPdtt = (Dv(P)*Du(ss) - Du(P)*Dv(ss))/det(J)

return (Dv(p)*Du(ss) - Du(p)*Dv(ss))/detJ; // dPdt
}

// polimorphic variations of the previous function
float DerivEx(float x; varying float ss, tt;)
{
float detJ = Du(ss)*Dv(tt) - Du(tt)*Dv(ss);
return (Du(x)*Dv(tt) - Dv(x)*Du(tt))/detJ;
}
vector DerivEx(vector vec; varying float ss, tt;)
{
float detJ = Du(ss)*Dv(tt) - Du(tt)*Dv(ss);
return (Du(vec)*Dv(tt) - Dv(vec)*Du(tt))/detJ;
}
color DerivEx(color c; varying float ss, tt;)
{
float detJ = Du(ss)*Dv(tt) - Du(tt)*Dv(ss);
return (Du(c)*Dv(tt) - Dv(c)*Du(tt))/detJ;
}


//// DETERMINANT OF THE JACOBIAN
//float det = Du(ss)*Dv(tt) - Du(tt)*Dv(ss);
//
//vector aniso_axis;
//
//// MULTIPLY BY THE INVERSE OF THE JACOBIAN
//if ( swap_axis > 0 )
//aniso_axis = normalize((Dv(P)*Du(ss) - Du(P)*Dv(ss))/det); // dPdt
//else
//aniso_axis = normalize((Du(P)*Dv(tt) - Dv(P)*Du(tt))/det); // dPds
//
//return aniso_axis;
float aw_turbulence( point p; uniform float frequency, octaves, lacunarity, gain)

{

	varying float sum = 0, amp =1, fade;
	
	float ss = s;
	float tt = t;
	varying point pp = frequency * transform("shader" , p);
	uniform float i;
	
	float filterWidthSS= abs(Du(ss)*du)+abs(Dv(ss)*dv);
	float filterWidthTT= abs(Du(tt)*du)+abs(Dv(tt)*dv);
	
	//float filterWidth = filterWidthSS>filterWidthTT?
	//filterWidthSS:filterWidthTT;
	float filterWidth = max(filterWidthSS, filterWidthTT);
	float val;
	
	for(i =0; i < octaves; i +=1) {
	
		fade =smoothstep(0.2, 0.3, filterWidth*gain);
		val = (1-fade) *(float noise(pp*gain)) + fade *0.5;
	
		sum += amp * ((val-0.5)/gain);
	amp *= gain; pp *= lacunarity;
	
	}
	
	return sum;

}

