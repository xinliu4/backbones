import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

min_z0_m200count = 2000
output_file_name = 'AlphaQ_formationtimes.txt'
particlemass = 1.14828e+9 # AlphaQ

subfiles_count = 8
total_hdf5_files = 256
file_edges = np.array([total_hdf5_files/subfiles_count*i for i in range(subfiles_count+1)])

# Mass definitions

# z_M(f_M): halo's main progenitor reaches a fraction F_M of the halo's final mass
F_M_threshold = np.array([0.01, 0.04, 0.05, 0.1, 0.25, 0.3, 0.5, 0.75, 0.8, 0.85, 0.9])

# z_T(M_T): main progenitor first reaches fixed mass
M_T_threshold = np.array([2.1e11, 1.*10**11.5, 1e12, 3.7e12, 1e14])/particlemass

# zW02 and zT04
formation_a_guess = 0.5
formation_p_guess = 0. # corresponds to W02, see T04 pg 5 right column

def zW02_form(ac, a0, M0, a):
	S = 2
	return M0*np.exp(-ac*S*(a0/a-1))

def zT04_form(ac, a0, M0, a, p):
	S = 2
	alpha = ac*S/a0
	return M0*(a/a0)**p * np.exp(alpha*(1-a0/a))

# z_vpeak(f_V): vpeak first reaches a fraction f_V of its maximum value
f_V_threshold = np.array([0.6, 0.8, 0.85, 0.9, 1])

# zvmax: vr200 reaches its maximum value
# zvvir: vr200 first reaches its z=0 value
# zZ03: maximum difference of circular velocity from expansion rate
def Z03_form(vh, z):
	gamma = -1./4.
	H_hubble = H_hubble_fcn(z)
	return np.log(vh) - gamma*np.log(H_hubble)

def H_hubble_fcn(redshift):
	params_list = np.array([0.220, 0.02258, 0.0, 0.71, -1.0, 0.0, 2.726, 0., 0.])
	temp_omega_cdm = params_list[0]
	temp_deut = params_list[1]
	temp_omega_nu = params_list[2]
	temp_hubble = params_list[3]   
	temp_w_de = params_list[4]
	temp_wa_de = params_list[5]
	temp_tcmb = params_list[6]
	temp_neff_massless = params_list[7]
	temp_neff_massive = params_list[8]
	H2 = hubble_sqrd(redshift, temp_omega_cdm, temp_deut/temp_hubble/temp_hubble, temp_omega_nu, temp_hubble, temp_tcmb, temp_neff_massless, temp_neff_massive, temp_w_de, temp_wa_de )
	return np.sqrt(H2)*100.*params_list[3]

def hubble_sqrd(redshift, Omega_cdm0, Omega_baryon0, Omega_nu0, h0, T_cmb0, N_eff_massless, N_eff_massive, w0, wa):
	m_Omega_cdm0 = Omega_cdm0
	m_Omega_baryon0 = Omega_baryon0
	m_Omega_nu0 = Omega_nu0
	m_h0 = h0
	m_T_cmb0 = T_cmb0
	m_N_eff_massless = N_eff_massless
	m_N_eff_massive = N_eff_massive
	m_w0 = w0
	m_wa = wa

	m_Omega_cb0 = m_Omega_cdm0 + m_Omega_baryon0
	m_Omega_matter0 = m_Omega_cb0 + m_Omega_nu0

	m_Omega_rad0 = 2.471e-5*(m_T_cmb0/2.725)**4.0/m_h0**2.0
	m_f_nu_massless = m_N_eff_massless*7.0/8.0*(4.0/11.0)**(4.0/3.0)
	m_f_nu_massive = m_N_eff_massive*7.0/8.0*(4.0/11.0)**(4.0/3.0)

	a = 1./(1.+redshift)
	mat = m_Omega_nu0/a**3.
	rad = m_f_nu_massive*m_Omega_rad0/a**4.
	Omega_nu_massive = (mat>=rad)*mat + (rad>mat)*rad
	return m_Omega_cb0/a**3. + (1.0 + m_f_nu_massless)*m_Omega_rad0/a**4. + Omega_nu_massive + (1.0 - m_Omega_matter0 - (1.0 + m_f_nu_massless)*m_Omega_rad0) *a**(-3.*(1.+m_w0+m_wa))*np.exp(-3.*m_wa*(1.-a))
 
#z_merge(Delta): the last time halo mass changed by a factor of at least (1+Delta) in one step
merge_Delta = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7])

#z_alpha Mstar(alpha): the last time halo mass crossed above 1000*Mstar(z)
Mstar_alpha_thresholds = [500., 1000., 1500.]
params_fit = [12.455, -2.015, 1.131, -6.278, 3.462, -13.266, 10.685]

def fit_fcn_klypin(z, A, B, b, C, c, D, d):
    y = z/(1.+z)
    return A + B*y**b + C*y**c + D*y**d

def Mstar_fcn(z):
    return 10**fit_fcn_klypin(z, params_fit[0], params_fit[1], params_fit[2], params_fit[3], params_fit[4], params_fit[5], params_fit[6])


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
	f.write('z_M(Ms/M200)\t')
	for M_T in M_T_threshold:
		f.write('z_T(%g)\t' %(M_T*particlemass))
	f.write('z_W02\tz_T04\t')
	for f_V in f_V_threshold:
		f.write('z_vpeak(%g)\t' %(f_V))
	f.write('z_vmax\tz_vvir\tz_Z03\t')
	for Delta in merge_Delta:
		f.write('z_merge(%g)\t' %Delta)
	for alpha in Mstar_alpha_thresholds:
		f.write('z_alphaMstar(%g)\t' %alpha)
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
				# M_{-2}(z=0)/M_{200}(z=0)
				if(z0_Ms < 0 or  min(backbone_mass >= z0_Ms)):
					z_M_MsM200 = -1
				elif (z0_Ms >= 0 and min(backbone_mass) < z0_Ms):
					 for node_i in range(1, len(backbone)):
						if backbone_mass[node_i-1] < z0_Ms and backbone_mass[node_i] >= z0_Ms:
							interp_mass = interp1d([backbone_mass[node_i-1], backbone_mass[node_i]],[backbone_redshift[node_i-1], backbone_redshift[node_i]])
							z_M_Ms_M200 = interp_mass(z0_Ms)
							break

				#print F_M_threshold
				#print z_M
				#print z0_Ms/z0_mass
				#print z_M_Ms_M200

				# fixed mass threshold
				z_T = np.ones_like(M_T_threshold)*-1
				for M_T_i in range(len(M_T_threshold)):
					M_T = M_T_threshold[M_T_i]
					if min(backbone_mass) < M_T:
						for node_i in range(1,len(backbone)):
							if backbone_mass[node_i-1] < M_T and backbone_mass[node_i] >= M_T:
								interp_mass = interp1d([backbone_mass[node_i-1], backbone_mass[node_i]],[backbone_redshift[node_i-1], backbone_redshift[node_i]])
								z_T[M_T_i] = interp_mass(M_T)
								break
					elif min(backbone_mass) >= M_T: z_T[M_T_i] = -1
			
				#print z_T	

				# fits to accretion history
				#	absolute_sigma only affects covariance matrix (i.e. cov values aren't one standard deviation, but parameters are unaffected)
				backbone_mass_err = np.sqrt(backbone_mass)
				zW02params, zW02cov = curve_fit(lambda x, ac: zW02_form(ac, 1, z0_mass, x), backbone_a, backbone_mass, p0=[formation_a_guess], sigma=backbone_mass_err) 
				zT04params, zT04cov = curve_fit(lambda x, ac, p: zT04_form(ac, 1, z0_mass, x, p), backbone_a, backbone_mass, p0=[formation_a_guess, formation_p_guess], sigma=backbone_mass_err)

				zW02 = 1./(zW02params[0])-1.
				zT04 = 1./(zT04params[0])-1.

				#print(backbone[:,1])
				#print(backbone_mass)
				#print(backbone_a)
				#print(backbone_redshift)
				#print(zW02)
				#print(zT04)
				#print(zT04params[1])

				# circular velocity definitions
				# zvpeak
				max_vpeak = np.max(backbone_vpeak)
				z_vpeak = np.ones_like(f_V_threshold)*-1
				for f_V_i in range(len(f_V_threshold)):
					f_V = f_V_threshold[f_V_i]
					if np.min(backbone_vpeak) < f_V*max_vpeak:
						for node_i in range(1,len(backbone)):
							if backbone_vpeak[node_i-1] < f_V*max_vpeak and backbone_vpeak[node_i] >= f_V*max_vpeak:
								interp_vpeak = interp1d([backbone_vpeak[node_i-1], backbone_vpeak[node_i]],[backbone_redshift[node_i-1], backbone_redshift[node_i]])
								z_vpeak[f_V_i] = interp_vpeak(f_V*max_vpeak)
								break
					elif np.min(backbone_vpeak) >= f_V*max_vpeak: z_vpeak[f_V_i] = -1

				# zvmax
				max_v200_ind = np.argmax(backbone_v200)
				zvmax = backbone_redshift[max_v200_ind]

				# zvvir
				z0_v200 = backbone_v200[-1]
				for node_i in range(1,len(backbone)):
					if backbone_v200[node_i-1] < z0_v200 and backbone_v200[node_i] >= z0_v200:
						interp_v200 = interp1d([backbone_v200[node_i-1], backbone_v200[node_i]],[backbone_redshift[node_i-1], backbone_redshift[node_i]])
						zvvir = interp_v200(z0_v200)
						break
				if(abs(zvvir) < 0.0001): zvvir = 0.	
	
				#print backbone_redshift
				#print backbone_v200
				#print z0_v200
				#print zvmax
				#print zvvir
	
				# zZ03
				z03_form_z = np.array([Z03_form(backbone_v200[i], backbone_redshift[i]) for i in range(len(backbone_redshift))])
				#print z03_form_z
				#print backbone_redshift
				#print H_hubble_fcn(backbone_redshift)
				#print backbone_v200
				zZ03 = backbone_redshift[np.argmax(z03_form_z)]

				# last major merger
				z_merge = np.ones_like(merge_Delta)*-1
				mass_change_fraction = np.array([backbone_mass[s_i]/backbone_mass[s_i-1] -1 for s_i in range(1, len(backbone_mass))])
				for Delta_i in range(len(merge_Delta)):
					Delta = merge_Delta[Delta_i]
					for node_i in reversed(range(len(backbone_mass)-1)):
						if(mass_change_fraction[node_i] >= Delta):
							z_merge[Delta_i] = backbone_redshift[node_i+1]
							break
				#print backbone_mass
				#print backbone_redshift
				#print mass_change_fraction
				#print merge_Delta
				#print z_merge

				# Mstar alpha
				z_Mstar_alpha = np.ones_like(Mstar_alpha_thresholds)*-1
				backbone_Mstar = Mstar_fcn(backbone_redshift)/particlemass
				for alpha_i in range(len(Mstar_alpha_thresholds)):
					alpha = Mstar_alpha_thresholds[alpha_i]
					if(backbone_mass[-1] >= alpha*backbone_Mstar[-1]): z_Mstar_alpha[alpha_i]=-2
					elif(backbone_mass[-1] < alpha*backbone_Mstar[-1]):
						for node_i in reversed(range(1,len(backbone_mass))):
							if(backbone_mass[node_i-1] >= alpha*backbone_Mstar[node_i-1] and backbone_mass[node_i] < alpha*backbone_Mstar[node_i]):
								interp_elem_1 = backbone_mass[node_i-1]-alpha*backbone_Mstar[node_i-1]
								interp_elem_2 =  backbone_mass[node_i]-alpha*backbone_Mstar[node_i]
								interp_massdiff=[interp_elem_2, interp_elem_1]
								interp_redshift = [backbone_redshift[node_i], backbone_redshift[node_i-1]]
								#print interp_massdiff
								#print interp_redshift

								interp_diff = interp1d(interp_massdiff, interp_redshift)
								z_Mstar_alpha[alpha_i] = interp_diff(0.)
								break
				#print backbone_Mstar
				#print backbone_mass
				#print backbone_redshift
				#print Mstar_alpha_thresholds
				#print z_Mstar_alpha	

				# OUTPUT
				f.write('%d\t%g\t%g\t' %(root_id, z0_mass, z0_concentration))
				for F_M_i in range(len(F_M_threshold)):
					f.write('%g\t' %z_M[F_M_i])
				f.write('%g\t' %z_M_Ms_M200)
				for M_T_i in range(len(M_T_threshold)):
					f.write('%g\t' %z_T[M_T_i])
				f.write('%g\t%g\t' %(zW02, zT04))
				for f_V_i in range(len(f_V_threshold)):
					f.write('%g\t' %z_vpeak[f_V_i])
				f.write('%g\t%g\t%g\t' %(zvmax, zvvir, zZ03))
				for Delta_i in range(len(merge_Delta)):
					f.write('%g\t' %z_merge[Delta_i])
				for alpha_i in range(len(Mstar_alpha_thresholds)):
					f.write('%g\t' %z_Mstar_alpha[alpha_i])
				f.write('\n')




