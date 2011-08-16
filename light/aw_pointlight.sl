
class aw_pointlight(
			   float intensity = 1;
			   color lightcolor = 1;
			   uniform float falloff = 1; // 2 = quadratic, 3 = cubic, quartic = 4
			   point from = point "shader" (0,0,0);
			   /* Shadow controls */
			   uniform color shadowColor = 0;
			   /* Transmission Controls */
			   uniform float doTransmission = 0;
			   uniform float transSamples = 16;
			   uniform float transMinSamples = 4;
			   uniform float transBlur = 0.01;
			   uniform float transBias = 0.01;
			   /* Extra utility params */
			   uniform float invertshadow = 0;
			   uniform float inverttrans = 0;
			   uniform float disableshadpass = 0;
			   output varying color _shadow;
               output varying color _cl_noshadow;
               output varying float _atten;
			   uniform float __nondiffuse = 0;
			   uniform float __nonspecular = 0;
			   string __category = "pxcLight,std,shad";

			   )
{
		
		public void light(output vector L; output color Cl;) {
				
		  /* Variables */
		  float dist;
          float atten = intensity;
		  color lc = 1;
		  vector Lrel;
          _cl_noshadow = 0;
		
		  illuminate( from ) {
			  if (N.L <= 0) {

			  	Lrel = vtransform("world", L);
			  	dist = length(L);

			  	atten *= 1 / (pow(dist, falloff)+1);
			  	_atten = atten;
		  		uniform float raydepth, s;
		  		rayinfo("depth", raydepth);
			  	s = (raydepth == 0) ? transSamples : 1;

			  	if (doTransmission != 0) {
			  	    _cl_noshadow += (atten * lightcolor);
			  	    _shadow = transmission(Ps, from, "samples", s,
			  	  		  "minsamples", transMinSamples,
			  	  		  "samplecone", transBlur,
			  	  		  "bias", transBias);
			  	    _shadow = inverttrans == 0 ? _shadow : 1 - lc;
			  	    Cl = mix(_cl_noshadow, atten*shadowColor, 1 - _shadow);
			  	} else {
			  	    _cl_noshadow += atten * lightcolor;
			  	    Cl = _cl_noshadow;
			  	}
			  }

		 } // illuminate
  }
} 
