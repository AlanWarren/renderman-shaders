/************************************************************************
 * aw_spotlight.sl - This is a standard spotlight with outputs for shadows
 * which aw_surface_std.sl imports and uses as an AOV
 *
 * Author: Alan Warren (bluemoonshine@gmail.com)
 *
 * $Revision: 1.1 $    $Date: 2009/03/12 $
 *
 ************************************************************************/

class aw_spotlight(
			   float intensity = 1;
			   color lightcolor = 1;
			   uniform float falloff = 1; // 2 = quadratic, 3 = cubic, quartic = 4
			   point from = point "shader" (0,0,0);
			   point to = point "shader" (0,0,1);
			   float  coneangle = radians(30);
			   float  conedeltaangle = radians(5);
			   float  beamDistribution = 2;
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
               output varying float _atten;
			   uniform float __nondiffuse = 0;
			   uniform float __nonspecular = 0;

			   )
{
		
		public void light(output vector L; output color Cl;) {
				
		  /* Variables */
		  float dist, cosangle;
          float atten = intensity;
		  color lc = 0;
          _cl_noshadow = 0;

          uniform vector A = (to - from) / length(to - from);
          uniform float cosoutside= cos(coneangle);
          uniform float cosinside = cos(coneangle-conedeltaangle);
		
		  illuminate(from,A,coneangle){
			dist = length(L);	  
			cosangle = L.A / dist;
			atten *= pow(cosangle, beamDistribution);
			atten *= smoothstep( cosoutside, cosinside, cosangle );
            _atten = atten;
            atten *= 1 / (pow(dist, falloff)+1);

			if (shadowMap != ""){
			  _shadow = shadow(shadowMap,Ps,"samples",shadowSamples,
							"blur",shadowBlur,
							"bias",shadowBias);
			  _shadow = invertshadow == 0 ? _shadow:1-_shadow;
              _cl_noshadow += (atten * lightcolor);
              lc = mix(lightcolor,shadowColor, _shadow);
			  Cl = atten  * lc;

			} else if (doTransmission != 0) {
              _cl_noshadow +=  (atten * lightcolor);
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
