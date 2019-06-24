import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

min_z0_m200count = 2000
output_file_name = 'AlphaQ_mass_fraction_formationtimes.txt'

subfiles_count = 8
total_hdf5_files = 256
file_edges = np.array([total_hdf5_files/subfiles_count*i for i in range(subfiles_count+1)])

# Mass definitions

# z_M(f_M): halo's main progenitor reaches a fraction F_M of the halo's final mass
F_M_threshold = np.linspace(0.01, 0.99, 99)
print F_M_threshold

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

with open(output_file_name, 'w') as f:
	f.write('rootID\tz0m200\tz0c200\t')
	for F_M in F_M_threshold:
		f.write('z_M(%g)\t' %F_M)
	f.write('\n')
	for backbone_file_i in range(subfiles_count):
		first_index = file_edges[backbone_file_i]
		last_index = file_edges[backbone_file_i+1]

		print("%d\t%d" %(first_index, last_index))

		backbones_all = np.genfromtxt("AlphaQ.%d_%d.haloTags_c_vpeak_vr200" %(first_index, last_index), skip_header=1)

		backbones = np.array([backbones_all[backbones_all[:,0]==root_id] for root_id in np.unique(backbones_all[:,0])])

		for backbone in backbones:
			if(backbone[-1][1] != 499):
				print ("wrong final step!")
				print backbone
				exit()
			if(min(backbone[:,0]) < 0 or min(backbone[:,1] < 0) or min(backbone[:,2])<0 or min(backbone[:,3])<0 or min(backbone[:,5])<0 or min(backbone[:,6]<0)):
				print ("error code where it shouldn't be!")
				exit()
			if(backbone[-1][3] >= min_z0_m200count):
				# CALCULATE FORMATION TIMES
				root_id = backbone[0][0]
				backbone_redshift = interp_z_step(backbone[:,1])
				backbone_a = interp_a_step(backbone[:,1])
				backbone_mass = backbone[:,3]
				z0_mass = backbone_mass[-1]
				z0_Ms = backbone[-1][7]
				backbone_concentration = backbone[:,4]
				z0_concentration = backbone_concentration[-1]
				backbone_vpeak = backbone[:,5]
				backbone_v200 = backbone[:,6]

				# fractions of final halo mass
				z_M = np.ones_like(F_M_threshold)*-1
				for F_M_i in range(len(F_M_threshold)):
					F_M = F_M_threshold[F_M_i]
					if min(backbone_mass) < F_M*z0_mass:
						for node_i in range(1,len(backbone)):
							if backbone_mass[node_i-1] < F_M*z0_mass and backbone_mass[node_i] >= F_M*z0_mass:
								interp_mass = interp1d([backbone_mass[node_i-1], backbone_mass[node_i]],[backbone_redshift[node_i-1], backbone_redshift[node_i]])
								z_M[F_M_i] = interp_mass(F_M*z0_mass)
								break
					elif min(backbone_mass) >= F_M*z0_mass: z_M[F_M_i] = -1

				# OUTPUT
				f.write('%d\t%g\t%g\t' %(root_id, z0_mass, z0_concentration))
				for F_M_i in range(len(F_M_threshold)):
					f.write('%g\t' %z_M[F_M_i])
				f.write('\n')




