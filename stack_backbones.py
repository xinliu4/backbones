import numpy as np
from scipy.interpolate import interp1d

finalmass = 5.e+13
bin_width = 0.1

tforms = [0.5, 1., 2., 3.]

filename_out = "backbones_binned_finalmass_%g.txt" %finalmass

#formation_time_index = 9 # for half-mass formation time
formation_time_index = 18 # for threshold mass of 2000 particles

min_z0_m200_count = 2000
particlemass = 1.14828e+9 # AlphaQ

step_z = dict(((499,0.), (487,0.024), (475,0.05), (464,0.075), (453,0.101), (442,0.128), (432,0.154), (421,0.184), (411,0.212),
               (401,0.242), (392,0.271), (382,0.304), (373,0.335), (365,0.364), (355,0.402), (347, 0.433646), (338,0.471), (331,0.502), (323,0.539), 
               (315,0.578), (307,0.618), (300,0.656), (293,0.695), (286,0.736), (279,0.779), (272,0.824), (266,0.865), (259, 0.91446), (253,0.959), 
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

subfiles_count = 8
total_hdf5_files = 256
file_edges = np.array([total_hdf5_files/subfiles_count*i for i in range(subfiles_count+1)])

formation_times = np.genfromtxt('AlphaQ_formationtimes.txt', skip_header=1)

hid_formation_dict = {}
for line in formation_times:
	hid_formation_dict[line[0]] = line[formation_time_index]

bb_mass = {}
bb_conc = {}
bb_count = {}
for tform in tforms:
	bb_mass[tform] = {}
	bb_conc[tform] = {}
	bb_count[tform] = {}
	for step in step_z:
		bb_mass[tform][step] = 0.
		bb_conc[tform][step] = 0.
		bb_count[tform][step] = 0

with open(filename_out, "w") as f:
	for backbone_file_i in range(subfiles_count):
		first_index = file_edges[backbone_file_i]
		last_index = file_edges[backbone_file_i+1]

		print("%d\t%d" %(first_index, last_index))

		backbones_all = np.genfromtxt("AlphaQ.%d_%d.haloTags_c_vpeak_vr200" %(first_index, last_index), skip_header=1)
		print "loaded"
		backbones = np.array([backbones_all[backbones_all[:,0]==root_id] for root_id in np.unique(backbones_all[:,0])])
		print "arrayed"

		for backbone in backbones:
			if(backbone[-1][1] != 499):
				print ("wrong final step!")
				print backbone
				exit()
			if(backbone[-1][3]*particlemass >= finalmass*(1.-bin_width) and backbone[-1][3]*particlemass <= finalmass*(1.+bin_width)):
				formation_time = hid_formation_dict[backbone[0][0]]
				for tform in tforms:
					if(formation_time >= tform*(1.-bin_width) and formation_time <= tform*(1.+bin_width)):
						for line in backbone:
							if(line[4] != -1):
								step = int(line[1])
								bb_mass[tform][step] += line[3]
								bb_conc[tform][step] += line[4]
								bb_count[tform][step] += 1
	
	for tform in tforms:	
		print tform
		for step in step_z:
			print step
			print bb_count[tform][step]
			if(bb_count[tform][step] > 0):					
				f.write("%g\t%g\t%g\t%g\t%g\n" %(tform, step, bb_mass[tform][step]/float(bb_count[tform][step])*particlemass, bb_conc[tform][step]/float(bb_count[tform][step]), bb_count[tform][step])) # formation bin, mass, c200, count			 	
