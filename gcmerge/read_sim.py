######################################################
## read simulation output and convert to FITS files ##
######################################################

import yt, glob, os
from . import inputs


def make_fits(file, property, outputdir, slice=False, fits=False, weight_field=None, image_res=None): 
    ds = yt.load(file)

    if "xray_emissivity" in property:
        emin = property.split('_')[1]
        emax = property.split('_')[2]
        xray_fields = yt.add_xray_emissivity_field(ds, emin, emax, table_type='apec', metallicity=0.3)

    if slice==True:
        p = yt.FITSSlice(ds, 'z', ("gas",property))
        p.writeto(outputdir+"/%s_slice_%d.fits" % (property, filenum), image_res=image_res)

    if fits==True:
        prj_fits = yt.FITSProjection(ds, "z", [("gas",property)], weight_field=weight_field, image_res=image_res)
        prj_fits.writeto(outputdir+"/%s_proj_%d.fits" %(property, files.index(file)))

    print( 'Done for snap', files.index(file))