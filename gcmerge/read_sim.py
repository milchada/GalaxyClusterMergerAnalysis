######################################################
## read simulation output and convert to FITS files ##
######################################################
import yt, glob, os, gc

def make_fits(files, filenum, pot=True, sb=True, temp=True):
    file = files[filenum]
    ds = yt.load(file)

    if sb or temp:
        xray_fields = yt.add_xray_emissivity_field(ds, 0.3, 7, table_type='apec', metallicity=0.3)

    if pot:
        p = yt.FITSSlice(ds, 'z', ("gas","gravitational_potential"))
        p.writeto("fitsfiles/potential/potential_slice_%d.fits" % (files.index(file)))
        print(" Potential done")
        del(p)
        gc.collect()

    if sb:
        prj_fits = yt.FITSProjection(ds, "z", ('gas','xray_photon_emissivity_0.3_7_keV'), weight_field='emission_measure')
        prj_fits.writeto("fitsfiles/photon_emissivity/xray_photon_emissivity_0.3_7_keV_proj_%d.fits" %(files.index(file)))
        print(" Xray SB done")
        del(prj_fits)
        gc.collect()

    if temp:
        prj_fits = yt.FITSProjection(ds, "z", ('gas','temperature'), weight_field='mazzotta_weighting')
        prj_fits.writeto("fitsfiles/temperature/temperature_proj_%d.fits" %(files.index(file)))
        print("Temp done")
        del(prj_fits)
        gc.collect()
        print( 'Done for snap', files.index(file))

