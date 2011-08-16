/************************************************************************
 * aw_mattAOV.sl - coshader for holdout matte's
 * you can call the aov with the following syntax:
 * color coshader(foo):matte
 *
 * Author: Alan Warren (bluemoonshine@gmail.com)
 *
 * $Revision: 1.1 $    $Date: 2011/07/20 $
 *
 ************************************************************************/

class aw_mattAOV(color input = 0;
              float mode = 0;
              string __category = "holdouts";
              output color matte;)
{
    public void out()
    {
        if (mode == 0)
            matte = input;
        else if (mode == 1)
            matte = color(0);
    }
}
