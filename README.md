########################## PLEASE NOTE: THIS BRANCH IS FOR MY PERSONAL USE ##########################

This package allows you to design, run, and analyse simulations of mergers of galaxy clusters.

1. Generate initial conditions

initial_conditions.py has everything you need to generate two halos with dark matter and non-radiative gas in hydrostatic equilibrium. You can specify the dark matter concentration, total mass, size and strength of gas cool core, and optionally add a BCG potential. 

The outputs are the particle files and profiles for each halo in HDF5 format. Use these to run GAMER, as documented in https://github.com/gamer-project/gamer/wiki

2. Make FITS files

Once the simulation has run, use make_fits.py to generate FITS files of the temperature (Mazzotta-weighted projection), photon emissivity (emission-weighted), projected density and a slice of the gravitational potential. By default, the projections are all done down the z axis, perpendicular to the merger plane, and the slice is along this plane. 

If you want off-axis images, specify the normal vector to your plane of interest in the make_fits function. 

#By default, these files are saved to fitsfiles/temperature, fitsfile/rhoproj, fitsfile/potential and fitsfile/photon_emissivity, so create these directories within the simulation directory before running make_fits.py#.

3. Find BCGs

bcg_dist.py finds the two potential minima in the gravitational potential slice, with a minimum separation as specified input. If ret_peaks = True, it returns the pixel coordinates of the two minima. Otherwise, it returns the distance between them in kpc. This is a sensitive measure of the dynamical time. 
*If available and well-identified*, use observed positions of BCGs to select best-fit time in the simulation. If not, see measuring X-ray observables below.

bcg_velocities.py selects DM particles within a specified radius of a potential minimum, and returns their average 3D velocity. 

4. Fit lensing profiles

concentration.py takes in a FITS image of projected density and a slice of gravitational potential in the plane of the merger. It identifies the two potential minima, and fits projected NFW halos around them. User may choose whether to fit cutouts around each center, or do a joint fit, i.e. fitting a sum of two NFWs. 

5. Measure velocity power spectrum

turbulence_Pk.py takes as input a simulation snapshot and coordinates for a region of interest. Within this region, it interpolates velocities onto a regular grid, and then computes the power spectrum of these velocities. Useful diagnostic of levels of bulk and turbulent movements.

6. Find and group edges

find_features.py uses Cauchy edge detection to identify sharp features in the images. It returns a list of pixel points associated with any edge in the image.
islands.py sets up object classes used to sort these pixels into contiguous features.
make_islands.py actually takes in the point list from find_features and sorts them into islands.

7. Fit arcs to features

For a given island (set of points), find_points_above.py selects a subset whose contrast with their background exceeds the input threshold. Thus, a high threshold means fewer points around the sharpest part of the feature, and a low threshold goes further out. 
fit_arcs.py sequentially lowers the threshold from 0.9 to 0.1 (normalised to the maximum value for the given set) and fits a circular arc to each resulting subset of points. It outputs the selected points, plus best fit radius, center of curvature, and reduced chi-sq of each fit. 

8. Measure profiles

make_profiles.py takes in two FITS image (temperature plus anything else) and an island generated in make_islands.py. The image center can be specified as either the X-ray center or the radius of curvature of the given feature. make_profile then makes a wedge from from the center to the feature and computes the radial profiles of two FITS images across it. Optionally, the input can include Matplotlib axes onto which the profiles are plotted.

### In Development ###

Unfortunately, the edge_detection algorithms find a few features that would not be deemed important by the human eye. I don't yet have a way to match the features to observed ones. However, the longest ones will be the ones observed. For example, the two shocks, the cold front, and the disrupted core of the main cluster will be the four longest edges. Currently, the best way to decide which is which is unfortunately manual inspection of the snapshots.
