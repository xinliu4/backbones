import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import binned_statistic
import math

formation_time = 'halfmass'
formation_time_index = 9
subfiles_count = 8
total_hdf5_files = 256
file_edges = np.array([total_hdf5_files/subfiles_count*i for i in range(subfiles_count+1)])

step_z = dict(((499,0.), (487,0.024), (475,0.05), (464,0.075), (453,0.101), (442,0.128), (432,0.154), (421,0.184), (411,0.212),
               (401,0.242), (392,0.271), (382,0.304), (373,0.335), (365,0.364), (355,0.402), (338,0.471), (331,0.502), (323,0.539), 
               (315,0.578), (307,0.618), (300,0.656), (293,0.695), (286,0.736), (279,0.779), (272,0.824), (266,0.865), (253,0.959), 
               (247,1.006), (241,1.055), (235,1.107), (230,1.152), (224,1.209), (219,1.258), (213,1.321), (208,1.376), (203,1.433),
               (198,1.494), (194,1.544), (189,1.610), (184,1.680), (180,1.738), (176,1.799), (171,1.880), (167,1.947), (163,2.020),
               (159,2.092), (155,2.170), (151,2.252), (148,2.317), (144,2.407), (141,2.478), (137,2.577), (134,2.655), (131,2.736),
               (127,2.851), (124,2.941), (121,3.040), (119,3.102), (116,3.205), (113,3.313), (110,3.427), (107,3.548), (105,3.631), 
               (102,3.763), (100,3.855), (97,4.000), (95,4.102), (92,4.262), (90,4.374), (88,4.492), (86,4.615), (84,4.743), 
               (81,4.95), (79,5.091), (77,5.242), (76,5.321), (74,5.484), (72,5.656), (70,5.837), (68,6.028), (67,6.128),
               (65,6.336), (63,6.556), (62,6.672), (60,6.913), (59,7.04), (57,7.306), (56,7.445), (54,7.739), (53,7.894), 
               (52,8.054), (50,8.393), (49,8.571), (48,8.757), (46,9.152), (45,9.361), (44,9.579), (43,9.806), (42,10.044) ))

step_z_chk = np.array([499,0.])
for i in range(500):
    try: step_z_chk = np.vstack((step_z_chk,[i,step_z[i]]))
    except KeyError: 
        j=0
interp_a_step = interp1d(step_z_chk[1::,0], 1./(1.+step_z_chk[1::,1]))
def interp_z_step(step):
    return 1./interp_a_step(step)-1.

formation_times_array = np.genfromtxt('AlphaQ_formationtimes.txt', skip_header=1)

formation_time_dictionary = {}
formation_mass_dictionary = {}
formation_concentration_dictionary = {}
for line in formation_times_array:
	formation_time_dictionary[line[0]] = line[formation_time_index] 

for backbone_file_i in range(subfiles_count):	
	first_index = file_edges[backbone_file_i]
	last_index = file_edges[backbone_file_i+1]

	print("%d\t%d" %(first_index, last_index))

	backbones_all = np.genfromtxt("AlphaQ.%d_%d.haloTags_c_vpeak_vr200" %(first_index, last_index), skip_header=1)
	print "backbones read"

	backbones = np.array([backbones_all[backbones_all[:,0]==root_id] for root_id in np.unique(backbones_all[:,0])])
	print "backbones arrayed"

	for backbone in backbones:
		if(backbone[-1][1] != 499):
			print ("wrong final step!")
			print backbone
			exit()
		backbone_rootID = backbone[-1][0]
		if(backbone_rootID in formation_time_dictionary):
			if(formation_time_dictionary[backbone_rootID] > 0):
				this_halo_tform = formation_time_dictionary[backbone_rootID]
				concentration_backbone = backbone[:,4]
				redshift_backbone = interp_z_step(backbone[:,1])
				#print redshift_backbone
				interpolate_formation_index = interp1d(redshift_backbone[::-1], [i for i in range(len(backbone))][::-1])
				#print this_halo_tform
				formation_index = interpolate_formation_index(this_halo_tform)
				#print formation_index
				concentration_below = concentration_backbone[math.floor(formation_index)]
				concentration_above = concentration_backbone[math.ceil(formation_index)]
				if(concentration_below != -1 and concentration_above != -1):
					interpolate_concentration = interp1d(redshift_backbone[::-1], concentration_backbone[::-1])
					formation_concentration_dictionary[backbone_rootID] = interpolate_concentration(this_halo_tform)
					mass_backbone = backbone[:,3]
					interpolate_mass = interp1d(redshift_backbone[::-1], mass_backbone[::-1])
					formation_mass_dictionary[backbone_rootID] = interpolate_mass(this_halo_tform)
	
formationtime_formationconc = np.array([[formation_time_dictionary[key], formation_concentration_dictionary[key], formation_mass_dictionary[key]] for key in formation_concentration_dictionary]) 

print formationtime_formationconc[:10]
			
#binned_time_conc = binned_statistic(formationtime_formationconc[:,0], formationtime_formationconc[:,1], 'mean', 30)[0]
#binned_time_conc_stdev = binned_statistic(formationtime_formationconc[:,0], formationtime_formationconc[:,1], np.std, 30)[0]
#binned_time = binned_statistic(formationtime_formationconc[:,0], formationtime_formationconc[:,0], 'mean', 30)[0]
#binned_counts = binned_statistic(formationtime_formationconc[:,0], formationtime_formationconc[:,1], 'count', 30)[0]
#binned_time_mass = binned_statistic(formationtime_formationconc[:,0], formationtime_formationconc[:,2], 'mean', 30)[0]
#binned_time_mass_stdev = binned_statistic(formationtime_formationconc[:,0], formationtime_formationconc[:,2], np.std, 30)[0]

binned_time_conc = binned_statistic(formationtime_formationconc[:,0], formationtime_formationconc[:,1], 'median', 30)[0]
binned_time_conc_25percent = binned_statistic(formationtime_formationconc[:,0], formationtime_formationconc[:,1], lambda y: np.percentile(y,25), 30)[0]
binned_time_conc_75percent = binned_statistic(formationtime_formationconc[:,0], formationtime_formationconc[:,1], lambda y: np.percentile(y,75), 30)[0]
binned_time = binned_statistic(formationtime_formationconc[:,0], formationtime_formationconc[:,0], 'median', 30)[0]
binned_counts = binned_statistic(formationtime_formationconc[:,0], formationtime_formationconc[:,1], 'count', 30)[0]
binned_time_mass = binned_statistic(formationtime_formationconc[:,0], formationtime_formationconc[:,2], 'median', 30)[0]
binned_time_mass_25percent = binned_statistic(formationtime_formationconc[:,0], formationtime_formationconc[:,2], lambda y:np.percentile(y,25), 30)[0]
binned_time_mass_75percent = binned_statistic(formationtime_formationconc[:,0], formationtime_formationconc[:,2], lambda y:np.percentile(y,75), 30)[0]

print binned_time_conc

with open("initial_concentration_"+formation_time+'_medians.txt', 'w') as f:
	for ii in range(len(binned_time)):
		f.write("%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n" %(binned_time[ii], binned_time_conc[ii], binned_time_conc_25percent[ii], binned_time_conc_75percent[ii], binned_time_mass[ii], binned_time_mass_25percent[ii], binned_time_mass_75percent[ii], binned_counts[ii]))
