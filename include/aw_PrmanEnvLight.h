/************************************************************************
 * aw_PrmanEnvlight.h - struct for AOV processing
 *
 * Author: Alan Warren (bluemoonshine@gmail.com)
 *
 * $Revision: 1.1 $    $Date: 2009/03/12 $
 *
 ************************************************************************/


#ifndef aw_PrmanEnvLight_h
#define aw_PrmanEnvLight_h

struct aw_PrmanEnvLight
{
    varying color refl_result = 0;
    varying color occ_result = 0;
    varying color ind_result = 0;
    varying color sss_result = 0;
    varying color spec_result = 0;
    varying color indspec_result = 0;
    varying color diff_result = 0;
    varying color diffdirect_result = 0;
    varying color diffind_result = 0;
    varying color color_result = 0;
    varying color gather_result = 0;
    varying color specshad_result = 0;
    varying color shadow_result = 0;
    varying color normal_result = 0;
    varying color fresKt_result = 0;
    varying color fresKr_result = 0;
    varying color glossy_result = 0;

    varying color sph_occ = 0;
    varying color sph_rad = 0;

    varying color reflAOV = 0;
    varying color occAOV = 0;
    varying color indAOV = 0;
    varying color sssAOV = 0;
    varying color diffAOV = 0;
    varying color diffdirAOV = 0;
    varying color indspecAOV = 0;
    varying color specAOV = 0;
    varying color colorAOV = 0;
    varying color gathAOV = 0;
    varying color shadowAOV = 0;
    varying color specshadAOV = 0;
    varying color normalAOV = 0;
    varying color fresKtAOV = 0;
    varying color fresKrAOV = 0;
    varying color diffindAOV = 0;
    varying color glossyAOV = 0;

    varying color res = 0;
    // items that get added
    // the rest are mult so
    // no special var required
    varying color indirect = 0;
    varying color sss = 0;

    public void sum()
    {
        reflAOV = refl_result;
        occAOV = occ_result;
        indAOV = ind_result;
        sssAOV = sss_result;
        diffAOV = diff_result;
        specAOV = spec_result;
        indspecAOV = indspec_result;
        diffdirAOV = diffdirect_result;
        gathAOV = gather_result;
        specshadAOV = specshad_result;
        shadowAOV = shadow_result;
        normalAOV = normal_result;
        colorAOV = color_result;
        fresKtAOV = fresKt_result;
        fresKrAOV = fresKr_result;
        diffindAOV = diffind_result;
        glossyAOV = glossy_result;
    }

}



#endif
