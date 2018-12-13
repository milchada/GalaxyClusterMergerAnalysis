#read gamer
import yt
import glob
import os

with open('../gamer_files.txt') as f:
	a=f.readlines()

sciencedir = a[1].split('\n')[0]
homedir = '/home/uchadaya'#'/home/fas/nagai/uc24'

ds = yt.load(sciencedir+'/Data_000001')
Msim = 6#e14Msun
Mobs = 4.15#e14Msun
scalingcoeff = Msim/Mobs

units_override = {"mass_unit":(scalingcoeff**1./3 * ds.mass_unit.value, ds.mass_unit.units),
		  "length_unit":(scalingcoeff**1./3 * ds.length_unit.value, ds.length_unit.units)}

def plot(sciencedir, property,zmin=None, zmax=None, startsnap=0, proj=False, fits=True, endsnap=None,weight_field=None):
        files=glob.glob(sciencedir+'/Data*')
        files.sort()
        print( len(files), ' snaps to go')
        if endsnap == None:
            endsnap = len(files)

        for file in files[startsnap:endsnap]:
                #wait how to scale x and y
           ds = yt.load(file,units_override=units_override)
	   if property == "xray_emissivity_0.3_7.0_keV":
		   xray_fields = yt.add_xray_emissivity_field(ds, 0.3, 7.0, table_type='apec', metallicity=0.3)
           if proj==True:
                   # ds.add_field('scaledprop', ds[property]*scaling
                   
                   p = yt.ProjectionPlot(ds, 'z', ("gas",property),width=(14,'Mpc'),
                        center='c',weight_field=weight_field)
                   p.set_zlim(property,zmin, zmax)
                   p.save()
           if fits==True:
                   prj_fits = yt.FITSProjection(ds, "z", [("gas",property)], weight_field=weight_field)
                   prj_fits.writeto("tempproj/%s_proj_%d.fits" %(property, files.index(file)))
           print( 'Done for snap', files.index(file))


#fantastic, this is working on local machine as is. all dependencies installed.
if __name__=='__main__':
	simid=2
	sciencedir = a[simid].split('\n')[0]
	plot(sciencedir, 'temperature',weight_field='mazzotta_weighting',proj=False,fits=True)

