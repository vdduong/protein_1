# backbone HN HA
import re
import math

class Point:
  def __init__(self, name = '', x=0.0, y=0.0, z=0.0):
    self.name = name
    self.x = x
    self.y = y
    self.z = z

dict_bb = dict()

file = open('/Users/James/protein_project/2KA0.bb','r')
for line in file:
	line = line.rstrip('\n\r')
	line = re.split(' +', line)
	#print line
	if line[0] == 'HA' or line[0] == 'H':
		name = line[0]+'_'+line[1]
		dict_bb[name] = Point(name, float(line[2]), float(line[3]), float(line[4]))
		#print name, dict_bb[name].x, dict_bb[name].y, dict_bb[name].z
file.close()

dict_distance_bb = dict()
dict_distance_bb_complet = dict()
for key_1 in dict_bb.keys():
	for key_2 in dict_bb.keys():
		distance = (dict_bb[key_1].x - dict_bb[key_2].x)**2 + \
			(dict_bb[key_1].y - dict_bb[key_2].y)**2 + \
			(dict_bb[key_1].z - dict_bb[key_2].z)**2
		distance = math.sqrt(distance)
		dict_distance_bb_complet[(key_1, key_2)] = distance
		if distance <= 5.5:
			dict_distance_bb[(key_1, key_2)] = distance
		else:
			dict_distance_bb[(key_1, key_2)] = float('inf')

def fw_algo(dict_distance):
	d = dict(dict_distance)
	for k in dict_bb.keys():
		for i in dict_bb.keys():
			for j in dict_bb.keys():
				d[(i,j)] = min(d[(i,j)], d[(i,k)] + d[(k,j)])
	return d

dict_distance_bb = fw_algo(dict_distance_bb)

for key in dict_distance_bb.keys():
	print key, dict_distance_bb[key], dict_distance_bb_complet[key]

distance_matrix_HH = dict_distance_bb

### probability chain for beta sheet structure

def prob_sheet_HN_HN_1(name_pair):
	global distance_matrix_HH
	mean_distance = 4.29
	tol_distance = 0.5
	return math.exp(-0.5*(distance_matrix_HH[name_pair]-mean_distance)**2/tol_distance**2)

def prob_sheet_HN_HN_2(name_pair):
	global distance_matrix_HH
	mean_distance = 6.64
	tol_distance = 0.61
	return math.exp(-0.5*(distance_matrix_HH[name_pair]-mean_distance)**2/tol_distance**2)

def prob_sheet_HA_i_HA_j(name_pair):
	'''probability of the distance between HA i and HA j in an anti parallel beta sheet structure'''
	global distance_matrix_HH
	mean_distance = 7.30
	tol_distance = 0.42
	return math.exp(-0.5*(distance_matrix_HH[name_pair]-mean_distance)**2/tol_distance**2)

### probability chain for helix alpha structure

def prob_helix_1(name_pair):
	mean_distance = 2.80
	tol_distance = 0.31
	return math.exp(-0.5*(distance_matrix_HH[name_pair]-mean_distance)**2/tol_distance**2)

def prob_helix_2(name_pair):
	mean_distance = 4.33
	tol_distance = 0.33
	return math.exp(-0.5*(distance_matrix_HH[name_pair]-mean_distance)**2/tol_distance**2)

def prob_helix_3(name_pair):
	mean_distance = 4.88
	tol_distance = 0.33
	return math.exp(-0.5*(distance_matrix_HH[name_pair]-mean_distance)**2/tol_distance**2)

def prob_helix_4(name_pair):
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

#############

def fragment_score(name_pair, relation):
	'''define fragment score based on network anchoring '''
	vertices_rest = set(dict_bb.keys())
	vertices_rest.remove(name_pair[0])
	vertices_rest.remove(name_pair[1])
	score = 0.0
	
	if relation == 'helix_1':
		for vertice in vertices_rest:
			score+=prob_helix_1((name_pair[0], vertice))*prob_helix_2((name_pair[1], vertice))
			score+=prob_helix_2((name_pair[0], vertice))*prob_helix_1((name_pair[1], vertice))
			score+=prob_helix_2((name_pair[0], vertice))*prob_helix_3((name_pair[1], vertice))
			score+=prob_helix_3((name_pair[0], vertice))*prob_helix_2((name_pair[1], vertice))
			score+=prob_helix_4((name_pair[0], vertice))*prob_helix_3((name_pair[1], vertice))
			score+=prob_helix_3((name_pair[0], vertice))*prob_helix_4((name_pair[1], vertice))
			score+=prob_helix_1((vertice, name_pair[0]))*prob_helix_2((vertice, name_pair[1]))
			score+=prob_helix_2((vertice, name_pair[0]))*prob_helix_1((vertice, name_pair[1]))
			score+=prob_helix_2((vertice, name_pair[0]))*prob_helix_3((vertice, name_pair[1]))
			score+=prob_helix_3((vertice, name_pair[0]))*prob_helix_2((vertice, name_pair[1]))
			score+=prob_helix_3((vertice, name_pair[0]))*prob_helix_4((vertice, name_pair[1]))
			score+=prob_helix_4((vertice, name_pair[0]))*prob_helix_3((vertice, name_pair[1]))
		return score
		
	elif relation == 'sheet_2':
		for vertice in vertices_rest:
			score+=prob_sheet_1((name_pair[0], vertice))*prob_sheet_1((vertice, name_pair[1]))
			score+=prob_sheet_1((name_pair[1], vertice))*prob_sheet_1((vertice, name_pair[0]))
		return score
	elif relation == 'sheet_1':
		for vertice in vertices_rest:
			score+=prob_sheet_1((name_pair[0], vertice))*prob_sheet_2((name_pair[1], vertice))
			score+=prob_sheet_2((name_pair[0], vertice))*prob_sheet_1((name_pair[1], vertice))
			score+=prob_sheet_1((vertice, name_pair[0]))*prob_sheet_2((vertice, name_pair[1]))
			score+=prob_sheet_2((vertice, name_pair[0]))*prob_sheet_1((vertice, name_pair[1]))
		return score
	elif relation == 'helix_2':
		for vertice in vertices_rest:
			score+=prob_helix_1((name_pair[0], vertice))*prob_helix_1((vertice, name_pair[1]))
			score+=prob_helix_1((vertice, name_pair[0]))*prob_helix_1((name_pair[1], vertice))
			
			score+=prob_helix_3((name_pair[0], vertice))*prob_helix_1((name_pair[1], vertice))
			score+=prob_helix_1((name_pair[0], vertice))*prob_helix_3((name_pair[1], vertice))
			score+=prob_helix_1((vertice, name_pair[0]))*prob_helix_3((vertice, name_pair[1]))
			score+=prob_helix_3((vertice, name_pair[0]))*prob_helix_1((vertice, name_pair[1]))			
			
			score+=prob_helix_4((vertice, name_pair[0]))*prob_helix_2((vertice, name_pair[1]))
			score+=prob_helix_2((vertice, name_pair[0]))*prob_helix_4((vertice, name_pair[1]))
			score+=prob_helix_2((name_pair[0], vertice))*prob_helix_4((name_pair[1], vertice))
			score+=prob_helix_4((name_pair[0], vertice))*prob_helix_2((name_pair[1], vertice))
		return score
	elif relation == 'helix_3':
		for vertice in vertices_rest:
			score+=prob_helix_1((name_pair[0], vertice))*prob_helix_2((vertice, name_pair[1]))
			score+=prob_helix_1((name_pair[1], vertice))*prob_helix_2((vertice, name_pair[0]))
			score+=prob_helix_2((name_pair[0], vertice))*prob_helix_1((vertice, name_pair[1]))
			score+=prob_helix_1((name_pair[1], vertice))*prob_helix_2((vertice, name_pair[0]))
			
			score+=prob_helix_4((name_pair[0], vertice))*prob_helix_1((name_pair[1], vertice))
			score+=prob_helix_4((name_pair[1], vertice))*prob_helix_1((name_pair[0], vertice))
			score+=prob_helix_4((vertice, name_pair[0]))*prob_helix_1((vertice, name_pair[1]))
			score+=prob_helix_4((vertice, name_pair[1]))*prob_helix_1((vertice, name_pair[0]))
		return score
	elif relation == 'helix_4':
		for vertice in vertices_rest:
			score+=prob_helix_1((name_pair[0], vertice))*prob_helix_3((vertice, name_pair[1]))
			score+=prob_helix_1((name_pair[1], vertice))*prob_helix_3((vertice, name_pair[0]))
			
			score+=prob_helix_2((name_pair[0], vertice))*prob_helix_2((vertice, name_pair[1]))
			
			score+=prob_helix_3((name_pair[0], vertice))*prob_helix_1((vertice, name_pair[1]))
			score+=prob_helix_3((name_pair[1], vertice))*prob_helix_1((vertice, name_pair[0]))
		return score
	

dict_distance_bb = fw_algo(dict_distance_bb)
class Root:
	def __init__(self, root_HN, root_HA):
		self.root_HN = root_HN
		self.root_HA = root_HA

nb = 0
list_nb = list()
set_root = set()
dict_root = dict()

# connection HN_HA intraresidue is located by 2.9 angstrom and tolerance of 0.25

for key in dict_distance_bb.keys():
	root_1 = re.split('_', key[0])[0]
	root_2 = re.split('_', key[1])[0]
	nb_2 = re.split('_', key[1])[1]
	if abs(dict_distance_bb[key] - 2.9) < 0.25 and root_1 == 'HN' and root_2 == 'HA': #and key[1] == 'HA_113':
		#print key, dict_distance_bb[key]
		dict_root[(key[0], key[1])] = Root(key[0], key[1])
		nb+=1

for key in dict_root.keys():
	print key, dict_root[key].root_HN, dict_root[key].root_HA

# given the weight of 2 to alpha_helix structure, 2 to beta_sheet structure and 1 to unstructured or loop

def relation_2_roots(root_1, root_2):
	pass
