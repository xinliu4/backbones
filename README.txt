This folder calculates halo formation times from the backbone mass accretion history. 

#### Backbones.cpp and Makefile ####
Combines sodproperties, sodpropertybins, and merger trees files to compute halo backbones. 

Output: a set of files with halo root ID (the halo ID of its z=0 descendant), step ID, halo ID at that step, halo mass (in counts) at that step, concentration at that step, peak circular velocity at that step, v200 at that step, and the mass enclosed in the scale radius at that step (where r_s is calculated from the fit concentration).

Parameters:
	file paths: one for hdf5 merger trees files, one for sodproperties files, and one for sodpropertybins files. Both are stored in tmp_string to loop over steps, but search for ".sodproperties", ".sodpropertybins", and ".hdf5" to modify the file paths. 
	total_hdf5_files is the number of hdf5 files in the merger tree files directory. The code will break these up into subfiles_count chunks.
	min_sod_count is the minimum halo size to record the mass and concentration of a halo in the haloproperties file. 
	min_z0_count_threshold is the minimum halo size (in number of particles) at z=0 for which backbones will be calculated. 2000 is good since smaller halos will not have concentrations.
	numbins is the number of radial bins in the sod profiles; this is 20 for all existing sod profiles.
	G_gravity is used to calculate circular velocities.	
	particlemass should be modified according to the simulation (1.14828e+9 for Alpha Quadrant, 1.05114e+8 for Q Continuum, and 1.848e+9 for Outer Rim)
	catalog_steps is the list of steps at which to add to the backbone. The count_steps=100 steps currently in this array are everything for Alpha Quadrant, Q Continuum, and Outer Rim, but Delta Quadrant will have more.

Command-line argument: an integer (from 0 to subfiles_count) to identify the subset of hdf5 merger trees files for which to calculate backbones. To do the entire catalog, you'll need to run this with every value up to subfiles_count, as done in run_script.sh.

#### run_script.sh ####
runs backbones.exe. 

Example command:
qsub -n 1 -t 60 ./run_script.sh



PYTHON
Once the backbones have been calculated, there are a few python files to calculate formation times and other quantities.

#### backbone_tforms.py ####
Calculates all formation times listed in section 4 except z_NFW.

Output: text file with columns for z=0 haloID, z=0 halo mass, z=0 halo concentration, and all formation times.

Parameters: 
	min_z0_m200count is the minimum halo size (in number of particles) at z=0 for which formation times will be calculated. 2000 is good since smaller halos will not have concentrations.
	output_file_name is the text file where formation times and final concentrations will be written
	particlemass should be modified according to the simulation (1.14828e+9 for Alpha Quadrant, 1.05114e+8 for Q Continuum, and 1.848e+9 for Outer Rim)
	subfiles_count and total_hdf5_files should be the same as in backbones.cpp.
	F_M_threshold is the array of fractions F_M for which formation times will be calculated for each backbone
	M_T_threshold is the array of threshold masses M_T for which formation times will be calculated for each backbone
	formation_a_guess and formation_p_guess are initial parameters for fitting zW02 and zT04 forms to formation histories; could be changed to test robustness of the fit
	f_V_threshold is the array of fractions f_V for which formation times z_vpeak will be calculated
	merge_Delta is the array of fractions Delta for which formation times z_merge will be calculated
	Mstar_alpha_thresholds is the array of values alpha for which formation times z_alpha_Mstar will be calculated. The two functions fit_fcn_klypin and Mstar_fcn give Mstar as a function of redshift for our cosmology (so need to be changed only if the cosmological parameters change).




#### mass_fraction_tforms.py ####
Similar to backbone_tforms.py, but only calculates formation times z_M(f_M) for different fractions f_M of the halo's final mass.
Used for Fig. 4 (as of March 19, 2019). 

Output: text file with columns for z=0 haloID, z=0 halo mass, z=0 halo concentration, and formation times z_M(f_M)

Parameters: 
	min_z0_m200count is the minimum halo size (in number of particles) at z=0 for which formation times will be calculated. 2000 is good since smaller halos will not have concentrations.
	output_file_name is the text file where formation times and final concentrations will be written
	subfiles_count and total_hdf5_files should be the same as in backbones.cpp.
	F_M_threshold is the array of fractions F_M for which formation times will be calculated




#### individual_backbones_select.py ####
This is an example of a way to select out an individual backbone--for example, like Fig. 6, to plot the mass and concentration evolution for a single halo as a function of redshift. Currently, it prints the backbones of the 10 smallest and 10 largest halos (based on z=0 mass) that formed before z=2.

Output: root ID (halo ID at redshift 0), step ID, halo mass at that step (in counts), halo concentration at that step

Parameters:
	formation_time_index selects which definition of formation time to use when deciding the halo's age (index is the column in the output of backbone_tforms.py)
	output_file_name is the text file where formation times and final concentrations will be written
	particlemass should be modified according to the simulation (1.14828e+9 for Alpha Quadrant, 1.05114e+8 for Q Continuum, and 1.848e+9 for Outer Rim)
	subfiles_count and total_hdf5_files should be the same as in backbones.cpp.
	formation_times imports the output of backbone_tforms.py




#### initial_concentration_mass_avg.py ####
Calculates the average concentration of halos forming at each redshift (e.g. Fig. 9).

Output: redshift, median concentration of halos forming at that redshift, 25th percentile of concentration, 75th percentile of concentration, median mass of halos forming at that redshift, 25th percentile of mass, 75th percentile of mass, number of halos forming at that redshift

Parameters:
	formation_time is a string for the type of formation time used. Only appears in the output file name.
	formation_time_index selects which definition of formation time to use when deciding the halo's age (index is the column in the output of backbone_tforms.py)
	subfiles_count and total_hdf5_files should be the same as in backbones.cpp.
	formation_times_array imports the output of backbone_tforms.py (change this file path if that output file path is modified)
	backbones_all imports the output of backbones.exe (change this file path if that output file path is modified)	




#### stack_backbones.py ####
An example of binning mass accretion histories by final halo mass and formation time (e.g. Figure 8). For a chosen final mass and set of formation redshifts, outputs the average mass and concentration histories along the average backbone.

Output: formation redshift, step along the average backbone, average mass at that step, average concentration at that step, number of halos identified at that step and used to compute the average

Parameters:
	finalmass is the z=0 mass bin; bin_width is the width of the bin. So the bin will be finalmass*(1.-bin_width) <= m200(z=0) <= finalmass*(1.+bin_width)
	tforms is the set of formation times for which to calculate stacked backbones. Width is again bin_width. 
	filename_out is the output file name.
	formation_time_index selects which definition of formation time to use when deciding the halo's age (index is the column in the output of backbone_tforms.py)
	min_z0_m200count is the minimum halo size (in number of particles) at z=0 for which formation times will be calculated. 2000 is good since smaller halos will not have concentrations.
	particlemass should be modified according to the simulation (1.14828e+9 for Alpha Quadrant, 1.05114e+8 for Q Continuum, and 1.848e+9 for Outer Rim)
	subfiles_count and total_hdf5_files should be the same as in backbones.cpp.
	formation_times imports the output of backbone_tforms.py
	backbones_all imports the output of backbones.exe (change this file path if that output file path is modified)	
