#######
# defining different distance probability for the alpha helix secondary structure 
def prob_helix_1(name_pair):
	'''probability of the distance between HN i and HN i+1'''
	mean_distance = 2.80
	tol_distance = 0.31
	return math.exp(-0.5*(distance_matrix_HH[name_pair]-mean_distance)**2/tol_distance**2)

def prob_helix_2(name_pair):
	'''probability of the distance between HN i and HN i+2'''
	mean_distance = 4.33
	tol_distance = 0.33
	return math.exp(-0.5*(distance_matrix_HH[name_pair]-mean_distance)**2/tol_distance**2)

def prob_helix_3(name_pair):
	'''probability of the distance between HN i and HN i+3'''
	mean_distance = 4.88
	tol_distance = 0.33
	return math.exp(-0.5*(distance_matrix_HH[name_pair]-mean_distance)**2/tol_distance**2)

def prob_helix_4(name_pair):
	'''probability of the distance between HN i and HN i+4'''
	mean_distance = 6.35
	tol_distance = 0.93
	return math.exp(-0.5*(distance_matrix_HH[name_pair]-mean_distance)**2/tol_distance**2)

### alpha_helix
def prob_HA_HA_i_2(name_pair):
	'''probability of the distance between H alpha i and i+2'''
	mean_distance = 6.61
	tol_distance = 0.26
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HA_HA_i_3(name_pair):
	'''probability of the distance between H alpha i and i+3'''
	mean_distance = 5.62
	tol_distance = 0.5
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HA_HA_i_4(name_pair):
	'''probability of the distance between H alpha i and i+4'''
	mean_distance = 6.68
	tol_distance = 1.15
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HA_HB_i_1(name_pair):
	'''probability of the distance between H alpha i and HB i+1'''
	mean_distance = 6.04
	tol_distance = 0.34
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HA_HB_i_2(name_pair):
	'''probability of the distance between H alpha i and HB i+2'''
	mean_distance = 6.09
	tol_distance = 0.63
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HA_HB_i_3(name_pair):
	'''probability of the distance between H alpha i and HB i+3'''
	mean_distance = 3.89
	tol_distance = 1.32
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HA_HB_i_4(name_pair):
	'''probability of the distance between H alpha i and HB i+4'''
	mean_distance = 6.26
	tol_distance = 1.29
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HA_HN_i_2(name_pair):
	'''probability of the distance between H alpha 1 and HN i+2'''
	mean_distance = 4.44
	tol_distance = 0.35
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HA_HN_i_3(name_pair):
	'''probability of the distance between H alpha 1 and HN i+3'''
	mean_distance = 3.55
	tol_distance = 0.43
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HA_HN_i_4(name_pair):
	'''probability of the distance between H alpha 1 and HN i+4'''
	mean_distance = 4.29
	tol_distance = 0.5
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HB_HA_i_1(name_pair):
	'''probability of the distance between HB i and HA i+1'''
	mean_distance = 4.71
	tol_distance = 0.37
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HB_HA_i_2(name_pair):
	'''probability of the distance between HB i and HA i+2'''
	mean_distance = 7.87
	tol_distance = 0.36
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HB_HA_i_3(name_pair):
	'''probability of the distance between HB i and HA i+3'''
	mean_distance = 7.94
	tol_distance = 0.65
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HB_HA_i_4(name_pair):
	'''probability of the distance between HB i and HA i+4'''
	mean_distance = 7.94
	tol_distance = 1.1
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HB_HB_i_1(name_pair):
	'''probability of the distance between HB i and HB i+1'''
	mean_distance = 5.72
	tol_distance = 0.60
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HB_HB_i_2(name_pair):
	'''probability of the distance between HB i and HB i+2'''
	mean_distance = 7.66
	tol_distance = 0.6
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HB_HB_i_3(name_pair):
	'''probability of the distance between HB i and HB i+3'''
	mean_distance = 6.16
	tol_distance = 1.11
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HB_HB_i_4(name_pair):
	'''probability of the distance between HB i and HB i+4'''
	mean_distance = 6.90
	tol_distance = 1.44
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HB_HN_i_1(name_pair):
	'''probability of the distance between HB i and HN i+1'''
	mean_distance = 3.36
	tol_distance = 0.4
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HB_HN_i_2(name_pair):
	'''probability of the distance between HB i and HN i+2'''
	mean_distance = 5.46
	tol_distance = 0.38
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HB_HN_i_3(name_pair):
	'''probability of the distance between HB i and HN i+3'''
	mean_distance = 5.54
	tol_distance = 0.37
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HB_HN_i_4(name_pair):
	'''probability of the distance between HB and HN i+4'''
	mean_distance = 5.68
	tol_distance = 0.4
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HN_HA_i_1(name_pair):
	'''probability of the distance between HN i and HA i+1'''
	mean_distance = 5.27
	tol_distance = 0.17
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HN_HA_i_2(name_pair):
	'''probability of the distance between HN i and HA i+2'''
	mean_distance = 7.04
	tol_distance = 0.25
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HN_HA_i_3(name_pair):
	'''probability of the distance between HN i and HA i+3'''
	mean_distance = 7.51
	tol_distance = 0.58
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HN_HA_i_4(name_pair):
	'''probability of the distance between HN i and HA i+4'''
	mean_distance = 8.77
	tol_distance = 0.85
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HN_HB_i_1(name_pair):
	'''probability of the distance between HN i and HB i+1'''
	mean_distance = 5.39
	tol_distance = 0.51
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HN_HB_i_2(name_pair):
	'''probability of the distance between HN i and HB i+2'''
	mean_distance = 6.12
	tol_distance = 0.58
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HN_HB_i_3(name_pair):
	'''probability of the distance between HN i and HB i+3'''
	mean_distance = 5.92
	tol_distance = 0.67
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

def prob_HN_HB_i_4(name_pair):
	'''probability of the distance between HN i and HB i+4'''
	mean_distance = 7.20
	tol_distance = 0.70
	return math.exp(-0.5*(distance_matrix[name_pair]-mean_distance)**2/tol_distance**2)

###############
## defining probabilities for different distances in antiparallel beta sheet structure
## residue backbones are typed by (HN, HA) distance <= 3.0 angstrom with tolerance of 0.25 angstrom
## output of the algorithm should be the probability of two roots to be adjacent,given the distance matrix

def prob_a_beta(root_1,root_2):
	
	def P_i_j(root_1,root_2):
		'''suppose that root_1,root_2 at position i and j'''
		P_i_j = exp(-0.5*(distance_matrix[(root_1.HA,root_2.HA)]-7.30)**2./0.42**2.)*\
				exp(-0.5*(distance_matrix[(root_1.HN,root_2.HN)]-3.04)**2./0.62**2.)*\
				exp(-0.5*(distance_matrix[(root_1.HA,root_2.HN)]-4.98)**2./0.63**2.)*\
				exp(-0.5*(distance_matrix[(root_1.HN,root_2.HA)]-4.83)**2./0.33**2.)
		return P_i_j
		
	def P_i_j_plus_1(root_1, root_2):
		'''suppose that root_1, root_2 at position i and j+1'''
		P_i_j_plus_1 = exp(-0.5*(distance_matrix[(root_1.HA, root_2.HN)]-7.79)**2./0.45**2.)*\
					exp(-0.5*(distance_matrix[(root_1.HA, root_2.HA)]-6.31)**2./0.64**2.)*\
					exp(-0.5*(distance_matrix[(root_1.HN, root_2.HA)]-3.71)**2./0.48**2.)*\
					exp(-0.5*(distance_matrix[(root_1.HN, root_2.HN)]-5.05)**2./0.28**2.)
		return P_i_j_plus_1
	
	def P_i_minus_1_j(root_1, root_2):
		'''suppose that root_1, root_2 at position i-1 and j'''
		P_i_minus_1_j = exp(-0.5*(distance_matrix[(root_1.HA, root_2.HA)]-5.65)**2./0.51**2.)*\
					exp(-0.5*(distance_matrix[(root_1.HA, root_2.HN)]-4.68)**2./0.62**2.)*\
					exp(-0.5*(distance_matrix[(root_1.HN, root_2.HN)]-7.00)**2./0.76**2.)*\
					exp(-0.5*(distance_matrix[(root_1.HN, root_2.HA)]-8.27)**2./0.66**2.)
		return P_i_minus_1_j
	
	def P_i_minus_1_j_plus_1(root_1, root_2):
		'''suppose that root_1, root_2 at position i-1 and j+1'''
		P_i_minus_1_j_plus_1 = exp(-0.5*(distance_matrix[(root_1.HA, root_2.HA)]-2.68)**2./0.67**2.)*\
					exp(-0.5*(distance_matrix[(root_1.HA, root_2.HN)]-5.02)**2./0.54**2.)*\
					exp(-0.5*(distance_matrix[(root_1.HN, root_2.HA)]-4.91)**2./0.58**2.)*\
					exp(-0.5*(distance_matrix[(root_1.HN, root_2.HN)]-7.55)**2./0.61**2.)
		return P_i_minus_1_j_plus_1
	
	def P_i_minus_1_j_plus_2(root_1, root_2):
		'''suppose that root_1, root_2 at position i-1 and j+2'''
		P_i_minus_1_j_plus_2 = exp(-0.5*(distance_matrix[(root_1.HA, root_2.HA)]-6.32)**2./0.65**2.)*\
					exp(-0.5*(distance_matrix[(root_1.HA, root_2.HN)]-3.76)**2./0.69**2.)*\
					exp(-0.5*(distance_matrix[(root_1.HN, root_2.HA)]-7.75)**2./0.66**2.)*\
					exp(-0.5*(distance_matrix[(root_1.HN, root_2.HN)]-5.07)**2./0.60**2.)
		return P_i_minus_1_j_plus_2
	
	P_adj = 0.0

	for key_root in dict_roots.keys():
		root_ = dict_roots[key_root]
		P_adj+=P_i_j_plus_1(root_1, root_)+P_i_minus_1_j_plus_1(root_, root_2)
		P_adj+=P_i_j_plus_1(root_2, root_)+P_i_minus_1_j_plus_1(root_, root_1)
		P_adj+=P_i_j(root_1, root_)+P_i_minus_1_j(root_2, root_)
		P_adj+=P_i_j(root_2, root_)+P_i_minus_1_j(root_1, root_)
	
	return P_adj
	
# for parallel beta sheet secondary structure
def prob_p_beta(root_1, root_2):
	def prob_i_j(root_1, root_2):
		P_i_j = exp()
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
