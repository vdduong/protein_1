# reconstructing the model from full and exact distance matrix.
dict_coord = dict()
import numpy as np
import re
import math
import itertools

class Point:
	def __init__(self, name='', x=0.0, y=0.0, z=0.0):
		self.name = name
		self.x = x
		self.y = y
		self.z = z

fm = open('/Users/James/protein_project/full_matrix_c','r')
dict_distance = dict()
for line in fm:
	line = line.rstrip('\n\r')
	line = re.split(' +', line)
	dict_distance[(line[0], line[1])] = float(line[2])


def next_point(atom_5, atom_1, atom_2, atom_3, atom_4):
	"""coordinate of the fifth point, find x5, y5, z5
	given the distances between the fifth point and four defined points"""
	global dict_coord, dict_distance
	x1, y1, z1 = dict_coord[atom_1].x, dict_coord[atom_1].y, dict_coord[atom_1].z
	x2, y2, z2 = dict_coord[atom_2].x, dict_coord[atom_2].y, dict_coord[atom_2].z
	x3, y3, z3 = dict_coord[atom_3].x, dict_coord[atom_3].y, dict_coord[atom_3].z
	x4, y4, z4 = dict_coord[atom_4].x, dict_coord[atom_4].y, dict_coord[atom_4].z

	d54 = dict_distance[(atom_5, atom_4)]
	d53 = dict_distance[(atom_5, atom_3)]
	d52 = dict_distance[(atom_5, atom_2)]
	d51 = dict_distance[(atom_5, atom_1)]
	
	
	a = np.array([[2.0*(x3-x4), 2.0*(y3-y4), 2.0*(z3-z4)], \
				 [2.0*(x2-x4), 2.0*(y2-y4), 2.0*(z2-z4)], \
				 [2.0*(x1-x4), 2.0*(y1-y4), 2.0*(z1-z4)]])
	
	b = np.array([d54**2.0 - d53**2.0 + x3**2.0 + y3**2.0 + z3**2.0 - x4**2.0 - y4**2.0 - z4**2.0, \
				d54**2.0 - d52**2.0 + x2**2.0 + y2**2.0 + z2**2.0 - x4**2.0 - y4**2.0 - z4**2.0, \
				d54**2.0 - d51**2.0 + x1**2.0 + y1**2.0 + z1**2.0 - x4**2.0 - y4**2.0 - z4**2.0])

	solve = np.linalg.solve(a,b)
	x5, y5, z5 = solve[0], solve[1], solve[2]
	d54_ = math.sqrt((x5-x4)**2.0 + (y5-y4)**2.0 + (z5-z4)**2.0)
	d53_ = math.sqrt((x5-x3)**2.0 + (y5-y3)**2.0 + (z5-z3)**2.0)
	d52_ = math.sqrt((x5-x2)**2.0 + (y5-y2)**2.0 + (z5-z2)**2.0)
	d51_ = math.sqrt((x5-x1)**2.0 + (y5-y1)**2.0 + (z5-z1)**2.0)
	

	return x5, y5, z5


# order of construction
#'HG22_30', 'HD13_30', 'HG23_30', 'HG13_30', 'HD1_32', \
list_local = ['HG21_30', 'HD11_30', 'HA_31', 'HG12_30', 'HB_30', 'HA_30', 'HD11_26', \
	'HA_27', 'HA_26', 'HG_26', 'HD22_26', 'HB2_29', 'HB2_26', \
	'HD12_26', 'HE1_32', 'HD23_26', 'H_32', 'H_30', 'HD3_31', 'HG3_31', \
	'H_27', 'H_26', 'HD13_26', 'H_28', 'HG3_27', 'HD12_30', 'HB3_29', \
	'HG23_22', 'HB3_25', 'HB3_26', 'HB3_27', 'HD21_26', 'HB3_32', 'HB3_31', \
	'HZ_32', 'HG13_23', 'HG11_22', 'HG21_22', 'HG21_23', 'HG3_28', 'HD22_25', \
	'HA_29', 'HA_28', 'HA_25', 'HA_24', 'HA_23', 'HA_22', 'HB2_31', 'HB2_32', \
	'HG_25', 'HD2_32', 'HG12_22', 'HG12_23', 'HB2_24', 'HE2_29', 'HA_32', 'HB2_28', \
	'HB2_27', 'HB2_25', 'HD2_31', 'HB2_20', 'HE3_27', 'HG2_28', 'HB_22', 'HB_23', \
	'HG2_21', 'HG2_27', 'HG2_24', 'HZ2_27', 'HD12_25', 'HD3_27', 'HD3_29', 'HE2_32', \
	'HD2_29', 'HD2_20', 'HD2_24', 'HD2_27', 'HE3_29', 'HG2_31', 'HD23_25', 'H_21', \
	'H_20', 'H_23', 'H_22', 'H_25', 'H_24', 'H_29', 'HD13_25', 'HG3_21', 'HG3_24',\
	 'HG3_29', 'HG2_29', 'HG22_23', 'HG22_22', 'HE2_20', 'HE2_24', 'HE2_27', 'HB3_28', \
	 'HG11_23', 'HG23_23', 'HZ3_24', 'HZ3_27', 'HB3_20', 'HB3_21', 'HB3_24', 'HB2_21', \
	 'HE3_24', 'HZ2_24', 'HD11_25', 'HZ2_29', 'HD1_20', 'HD3_24', 'HA_21', 'HA_20', \
	 'HG13_22', 'HD21_25', 'HZ1_27', 'HZ1_24', 'HZ1_29', 'HZ3_29', 'HE1_20', 'HH_20']

dict_coord['HG22_30'] = Point('HG22_30', 12.082, 0.486, -2.391)
dict_coord['HD13_30'] = Point('HD13_30', 9.280, -1.935, -3.093)
dict_coord['HG23_30'] = Point('HG23_30', 13.106, 1.428, -3.450)
dict_coord['HG13_30'] = Point('HG13_30', 10.223, -1.249, -5.203)
dict_coord['HD1_32'] = Point('HD1_32', 9.037, 1.509, -6.286)

tol_dis = 0.5

for item in list_local:
#if 1:
	dict_cluster = dict()
	nb_cluster = 0
	for subset in itertools.combinations(dict_coord.keys(), 4):
		atom_1, atom_2, atom_3, atom_4 = subset[0], subset[1], subset[2], subset[3]
		x_next, y_next, z_next = next_point(list_local[0], atom_1, atom_2, atom_3, atom_4)
		if nb_cluster == 0:
			nb_cluster+=1
			dict_cluster[nb_cluster] = [(x_next, y_next, z_next)]
		else:
			for key_cluster in dict_cluster.keys():
				for item_cluster in dict_cluster[key_cluster]:
					distance_ = (item_cluster[0] - x_next)**2 + (item_cluster[1] - y_next)**2 + \
						(item_cluster[2] - z_next)**2
					if distance_ <= 0.25:
						dict_cluster[key_cluster].append((x_next, y_next, z_next))
						break
	for key_cluster in dict_cluster.keys():
		sum_x, sum_y, sum_z = 0.0, 0.0, 0.0
		for item_cluster in dict_cluster[key_cluster]:
			sum_x += item_cluster[0]
			sum_y += item_cluster[1]
			sum_z += item_cluster[2]
		count_ = len(dict_cluster[key_cluster])
		x_, y_, z_ = sum_x/count_, sum_y/count_, sum_z/count_
		print 'cluster nb %i :'%key_cluster, x_, y_, z_
	if len(dict_cluster.keys()) == 1:
		dict_coord[item] = Point(item, x_, y_, z_)
	else:
		pass
	
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
