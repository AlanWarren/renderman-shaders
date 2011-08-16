#define ALPHA 7

class aw_distantlight(
			   float intensity = 1;
			   color lightcolor = 1;
			   uniform float falloff = 1; // 2 = quadratic, 3 = cubic, quartic = 4
			   point from = point "shader" (0,0,0);
			   point to = point "shader" (0,0,1);
			   /* Shadow controls */
			   uniform color shadowColor = 0;
			   uniform string shadowMap    = "";
			   uniform float shadowBlur = 0.01; 
               uniform float shadowSamples = 32;
			   uniform float shadowBias = 0.01;
			   uniform float shadowWidth = 1;
			   uniform string shadowFilter = "gaussian";
			   /* Transmission Controls */
			   uniform float doTransmission = 0;
               uniform float shadowDist = 1000;
			   uniform float transSamples = 16;
			   uniform float transMinSamples = 4;
			   uniform float transBlur = 0.01;
			   uniform float transBias = 0.01;
			   /* Extra utility params */
			   uniform float invertshadow = 0;
			   uniform float inverttrans = 0;
			   string __category = "pxcLight,std,shad";
			   uniform float disableshadpass = 0;
			   output varying color _shadow;
               output varying color _cl_noshadow;
			   uniform float __nondiffuse = 0;
			   uniform float __nonspecular = 0;

			   )
{
		
		public void light(output vector L; output color Cl;) {
				
		  /* Variables */
		  float dist;
          float atten = intensity;
          vector A = to - from;
		  color lc = 0;
          _cl_noshadow = 0;

		
		  solar(A, 0.0) {
			dist = length(L);	  
            point PL = Ps + (shadowDist*normalize(-A));
			//cosangle = L.A / dist;

			if (shadowMap != ""){
			  _shadow = shadow(shadowMap,Ps,"samples",shadowSamples,
                            "width",shadowWidth,
							"blur",shadowBlur,
							"bias",shadowBias,
                            "filter", shadowFilter);
			  _shadow = invertshadow == 0 ? _shadow:1-_shadow;
              _cl_noshadow += lightcolor;
              Cl = mix(lightcolor,shadowColor, _shadow);
			  //Cl = atten  * lc;
			} else if (doTransmission != 0) {
              _cl_noshadow +=  lightcolor;
		      _shadow = transmission(Ps, from,"samples", transSamples,
		       				"minsamples", transMinSamples,
		       				"samplecone", transBlur,
		       				"bias", transBias);
		      _shadow = inverttrans == 0 ? _shadow: 1-lc;
              Cl = mix(_cl_noshadow, atten*shadowColor, 1-_shadow);
	
			}
            else {
			//Cl = atten * lightcolor; 
            _cl_noshadow += atten * lightcolor;
            Cl = _cl_noshadow;
            }
		 } // illuminate
  }
} 



