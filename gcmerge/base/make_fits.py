######################################################
## read simulation output and convert to FITS files ##
######################################################
import yt, glob, os, gc
import numpy as np 

def total_dens(field, data):
    #create a total density field
    return data['gas', 'density'] + data['particle_density_on_grid'] 

def make_fits(file, pot=True, sb=True, temp=True, rho=True):
    #read in simulation file
    ds = yt.load(file)
    #read in time; we use this to label outputs
    filenum = int(current_time.in_units("Gyr")*100)

    if sb or temp:
        #xray emissivity and temperature are not included by default
        #compute them here 
        xray_fields = yt.add_xray_emissivity_field(ds, 0.3, 7, table_type='apec', metallicity=0.3)

    if pot:
        #make sure not to remake any files
        if not glob.glob("fitsfiles/potential/potential_slice_%d.fits" % filenum):
            p = yt.FITSSlice(ds, 'z', ("gas","gravitational_potential"))
            p.writeto("fitsfiles/potential/potential_slice_%d.fits" % filenum)
            print(" Potential done")
            del(p)
            gc.collect()

    if sb:
        if not glob.glob("fitsfiles/photon_emissivity/xray_photon_emissivity_0.3_7_keV_proj_%d.fits" % filenum):
            prj_fits = yt.FITSProjection(ds, "z", ('gas','xray_photon_emissivity_0.3_7_keV'), weight_field=None)
            prj_fits.writeto("fitsfiles/photon_emissivity/xray_photon_emissivity_0.3_7_keV_proj_%d.fits" % filenum)
            print(" Xray SB done")
            del(prj_fits)
            gc.collect()

    if temp:
        if not glob.glob("fitsfiles/temperature/temperature_proj_%d.fits" % filenum):
            prj_fits = yt.FITSProjection(ds, "z", ('gas','temperature'), weight_field='mazzotta_weighting')
            prj_fits.writeto("fitsfiles/temperature/temperature_proj_%d.fits" % filenum)
            print("Temp done")
            del(prj_fits)
            gc.collect()

    if rho:
        if not glob.glob("fitsfiles/rhoproj/rhoproj_%d.fits" % filenum):
            ds.add_field(("gamer", "total_density"), 
            units="g/cm**3", function=total_dens, 
            sampling_type="cell")      

            prj_fits = yt.FITSProjection(ds, "z", 'total_density', weight_field=None)
            prj_fits.writeto("fitsfiles/rhoproj/rhoproj_%d.fits" % filenum)
            print("Dens done")
            del(prj_fits)
            gc.collect()

    print( 'Done for snap', filenum)

if __name__ == "__main__":
    files = glob.glob('Data*')
    files.sort()
    for file in files:
        make_fits(file)
