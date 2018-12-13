#read gamer
import yt
import glob
import os
from __init__ import *

with open(inputs['filenamedir']) as f:
	a=f.readlines()

sciencedir = a[int(inputs['simid'])].split('\n')[0]

ds = yt.load(inputs['sciencedir']+'/'+inputs['testdatafile'])

if rescale:
  # scalingcoeff = Msim/Mobs
  units_override = {"mass_unit":(rescale**1./3 * ds.mass_unit.value, ds.mass_unit.units),
		  "length_unit":(rescale**1./3 * ds.length_unit.value, ds.length_unit.units)}

def make_fits(sciencedir, property,zmin=None, zmax=None, startsnap=0, proj=False, fits=True, endsnap=None,weight_field=None):
        files=glob.glob(sciencedir+'/Data*')
        files.sort()
        print( len(files), ' snaps to go')
        if endsnap == None:
            endsnap = len(files)

        for file in files[startsnap:endsnap]:
          if rescale:
            ds = yt.load(file,units_override=units_override)
          else:
            ds = yt.load(file)

	   if property == "xray_emissivity_0.3_7.0_keV":
		   xray_fields = yt.add_xray_emissivity_field(ds, 0.3, 7.0, table_type='apec', metallicity=0.3)
           if proj==True:
            p = yt.ProjectionPlot(ds, 'z', ("gas",property),width=(14,'Mpc'),
                center='c',weight_field=weight_field)
            p.set_zlim(property,zmin, zmax)
            p.save()
           if fits==True:
            prj_fits = yt.FITSProjection(ds, "z", [("gas",property)], weight_field=weight_field)
            prj_fits.writeto("tempproj/%s_proj_%d.fits" %(property, files.index(file)))
           print( 'Done for snap', files.index(file))