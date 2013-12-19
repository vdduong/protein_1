import re
import math
import numpy as np
import random
import itertools



def fit(P,Q):
	'''	varies the distance between P and Q,
	and optimizes rotation for each step until a minimum is found	'''
	step_size = P.max(0)
	threshold = step_size*1e-5
	rmsd_best = kabsh(P,Q)
	while True:
		for i in range(3):
			temp = np.zeros(3)
			temp[i] = step_size[i]
			rmsd_new = kabsh(P + temp, Q)
			if rmsd_new < rmsd_best:
				rmsd_best = rmsd_new
				P[:,i]+=step_size[i]
			else:
				rmsd_new = kabsh(P-temp, Q)
				if rmsd_new < rmsd_best:
					rmsd_best = rmsd_new
					P[:,i]-=step_size[i]
				else:
					step_size[i]/=2
			if (step_size < threshold).all():
				break
	return rmsd_best
	
	
def kabsh(P,Q):
	""" Kabsh algorithm to align two point clouds """
	#print P
	C = np.dot(np.transpose(P), Q)
	V, S, W = np.linalg.svd(C)
	d = (np.linalg.det(V)*np.linalg.det(W)) < 0.0
	if (d):
		S[-1] = - S[-1]
		V[:,-1] = - V[:,-1]
	U = np.dot(V, W)
	P = np.dot(P, U)
	
	return rmsd(P, Q), P
	
def rmsd(V, W):
	''' return the rmsd from two sets of vectors V and W	'''
	D = len(V[0])
	N = len(V)
	rmsd = 0.0
	for v, w in zip(V, W):
		rmsd+=sum((v[i]-w[i])**2.0 for i in range(D)) 
	return np.sqrt(rmsd/N)
	
def centroid(X):
	'''	return the centroid from a vectorset X	'''
	C = sum(X)/len(X)
	return C
	
	
class Point:
	def __init__(self, name ='', x=0.0, y=0.0, z=0.0):
		self.name = name
		self.x = x
		self.y = y
		self.z = z

# dict_point defines the pdb coordinate of this protein fragment
dict_point = {}
file = open("/Users/James/protein_project/2KA0.pdb.1","r")
for line in file:
	line = line.rstrip('\n\r')
	line = re.split(" +", line)
	name = line[2] + '_'+line[5]
	if 'H' in name :#and 'O' not in name:
		coord_x = float(line[6])
		coord_y = float(line[7])
		coord_z = float(line[8])
		point = Point(name, coord_x, coord_y, coord_z)
		dict_point[name] = point
	else:
		continue


fm = open('/Users/James/protein_project/full_matrix_c','r')
g = dict()
g1 = dict()
for line in fm:
	line = line.rstrip('\n\r')
	line = re.split(' +', line)
	g1[(line[0], line[1])] = float(line[2])
	if float(line[2]) <= 5.5:
		g[(line[0], line[1])] = float(line[2])
	else:
		if line[0] == line[1]:
			g[(line[0], line[1])] = 0.0
		else:
			g[(line[0], line[1])] = float('inf')
			
			

vertices = set()
for item in g.keys():
	vertices.add(item[0])
	vertices.add(item[1])

# generate the adjacency matrix/ distance matrix from an incomplete one
# using Floyd Warshall all pairs of closest distance

def adj(vertices, dict_distance):
	"""
	converts the dictionary of distances to an adjacency matrix
	"""
	adj = {}
	for i in vertices:
		adj[i] = {}
		for j in vertices:
			try:
				adj[i][j] = dict_distance[(i,j)]
			except KeyError:
				# the distance from a node to itself is zero
				if i == j:
					adj[i][j] = 0.0
				# the distance from a node to an unconnected one is infinity
				else:
					adj[i][j] = float('inf')

def fw(vertices, dict_distance):
	"""
	run the FW algorithm
	"""
	d = dict(dict_distance)
	for k in vertices:
		for i in vertices:
			for j in vertices:
				d[(i,j)] = min(d[(i,j)], d[(i,k)] + d[(k,j)])
	return d

d = fw(vertices, g)
#d1 = fw(vertices, g)
d1 = g1
#for key in d.keys():
#	print key, d[key], g[key]

up_lim = dict()
low_lim = dict()
for key in d.keys():
	if d[key] <= 5.5:
		up_lim[key] = d[key]
		low_lim[key] = d[key]
	else:
		up_lim[key] = d[key]
		low_lim[key] = 1.3 # vdw radius
		
def triangle_ineq(vertices, dict_distance):
	vertices = list(vertices)
	for k in range(len(vertices)):
		for i in range(len(vertices)-1):
			for j in range(i+1, len(vertices)):
				if up_lim[(vertices[i], vertices[j])] > up_lim[(vertices[i], vertices[k])] + \
					up_lim[(vertices[k], vertices[j])]:
					up_lim[(vertices[i], vertices[j])] = up_lim[(vertices[i], vertices[k])] + \
					up_lim[(vertices[k], vertices[j])]
				if low_lim[(vertices[i], vertices[j])] < low_lim[(vertices[i], vertices[k])] - \
					up_lim[(vertices[k], vertices[j])]:
					low_lim[(vertices[i], vertices[j])] = low_lim[(vertices[i], vertices[k])] - \
					up_lim[(vertices[k], vertices[j])]
				elif low_lim[(vertices[i], vertices[j])] <= low_lim[(vertices[j], vertices[k])] - \
					up_lim[(vertices[k], vertices[j])]:
					low_lim[(vertices[i], vertices[j])] = low_lim[(vertices[j], vertices[k])] - \
					up_lim[(vertices[k], vertices[i])]

# applying the triangle inequality to the distance matrix
triangle_ineq(vertices, d)

#for key in d.keys():
#	print key, d[key], up_lim[key], low_lim[key], '%.4f'%((0.5*up_lim[key]**4.0 + 0.5*low_lim[key]**4.0)**0.25), \
#			g1[key]

dict_coord = {}

def constructing_4_points(list_4_points, dict_distance):
	""""
	constructs the first 4 points
	"""
	dict_coord = dict()
	atom_1 = list_4_points[0]
	x1, y1, z1 = 0.0, 0.0, 0.0
	dict_coord[atom_1] = Point(atom_1, x1, y1, z1)
	
	atom_2 = list_4_points[1]
	x2 = float('%.4f'%dict_distance[(atom_1, atom_2)])
	y2 = 0.0
	z2 = 0.0
	dict_coord[atom_2] = Point(atom_2, x2, y2, z2)
	
	atom_3 = list_4_points[2]
	d31 = dict_distance[(atom_1, atom_3)]
	d32 = dict_distance[(atom_2, atom_3)]
	x3 = float('%.4f'%((d31**2.0 - d32**2.0)/(2.0*x2) + x2/2.0))
	y3 = float('%.4f'%((d31**2.0 - x3**2.0)**(0.5)))
	z3 = 0.0
	dict_coord[atom_3] = Point(atom_3, x3, y3, z3)
	
	atom_4 = list_4_points[3]
	d41 = dict_distance[(atom_1, atom_4)]
	d42 = dict_distance[(atom_2, atom_4)]
	d43 = dict_distance[(atom_3, atom_4)]
	x4 = float('%.4f'%((d41**2.0 - d42**2.0)/(2.0*x2) + x2/2.0))
	y4 = float('%.4f'%((d42**2.0 - d43**2.0 - (x4-x2)**2.0 + (x4-x3)**2.0)/(2.0*y3) + y3/2.0))
	z4 = float('%.4f'%((d41**2.0 - x4**2.0 - y4**2.0)**0.5))
	dict_coord[atom_4] = Point(atom_4, x4, y4, z4)
	
	return dict_coord


def fifth_point(atom_5, list_4_points, dict_distance, dict_coord):
	"""coordinate of the fifth point, find x5, y5, z5
	given the distances between the fifth point and four defined points"""
	#global dict_coord
	atom_1, atom_2, atom_3, atom_4 = list_4_points[0], list_4_points[1], list_4_points[2], list_4_points[3]
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
	
	x5 = float('%.4f'%x5)
	y5 = float('%.4f'%y5)
	z5 = float('%.4f'%z5)

	return x5, y5, z5

D = [ ]
list_vertices = list(vertices)
list_c = list(vertices)
for key in list_vertices:
	x_d, y_d, z_d = dict_point[key].x, dict_point[key].y, dict_point[key].z
	coord_d = [x_d, y_d, z_d]
	D.append(np.array(coord_d))
D = np.array(D)
Dc = centroid(D)
D-=Dc


def iteration_procedure(d, number_iteration):
	global vertices, D
	list_c = list(vertices)
	list_vertices = list(vertices)
	dict_nb_model = {}
	dict_model = {}
	
	nb_model = 0
	nb_try = 0
	while nb_try <= number_iteration:
		while 1:
			try: 
				random.shuffle(list_c)
				list_4_points = list_c[:4]
				dict_coord_1 = constructing_4_points(list_4_points, d)
				break
			except ValueError as err:
				pass
			except ZeroDivisionError as err:
				pass
		
		for i in range(4, len(list_c)):
			next_point = list_c[i]
			x_n, y_n, z_n = fifth_point(next_point, list_4_points, d, dict_coord_1)
			dict_coord_1[next_point] = Point(next_point, x_n, y_n, z_n)
		
		V1 = []
		for key in list_vertices:
			x_1, y_1, z_1 = dict_coord_1[key].x, dict_coord_1[key].y, dict_coord_1[key].z
			coord_1 = [x_1, y_1, z_1]
			V1.append(np.array(coord_1))
		V1 = np.array(V1)
		V1c = centroid(V1)
		V1-=V1c
		
		rmsd_kabsh_local, V1 = kabsh(V1, D)
		dict_model[nb_try] = V1
		
		#if len(dict_model.keys()) == 0:
		#	nb_model +=1
		#	dict_nb_model[nb_model] = 1
		#	dict_model[nb_model] = V1
		#else:
		#	nb_found = 0
		#	for key in dict_model.keys():
		#		model_average = dict_model[key]
		#		rmsd_kabsh_local, V1 = kabsh(V1, model_average)
		#		if rmsd_kabsh_local <= 4.0:
		#			model_average = model_average*(dict_nb_model[nb_model]) + V1
		#			model_average = model_average/(dict_nb_model[nb_model] + 1)
		#			dict_nb_model[nb_model] += 1
		#			dict_model[nb_model] = model_average
		#			nb_found +=1
		#	if nb_found < 1:
		#		nb_model+=1
		#		dict_nb_model[nb_model] = 1
		#		dict_model[nb_model] = V1
		nb_try += 1
		
	#print 'there are in total %i models constructed'%nb_model
	return dict_model	
		


def corr(dict_model, d, number_iteration):
	global vertices, D
	dict_corr = dict()
	dict_count = dict()
	for key in d.keys():
		dict_corr[key] = set()
	for key in list_vertices:
		dict_count[key] = 0
		
	for i in range(0, len(dict_model.keys())):
		configuration = dict_model[i]
		dict_position = dict()
		for key in range(0, len(list_vertices)):
			dict_position[key] = 0
		
		for key_ in dict_model.keys():
			if key_ != i:
				other_configuration = dict_model[key_]
				rmsd_local, other_configuration = kabsh(other_configuration, configuration)
				for i in range(0, D.shape[0]):
					if (configuration[i,:][0] - other_configuration[i,:][0])**2 +\
					(configuration[i,:][1] - other_configuration[i,:][1])**2 +\
					(configuration[i,:][2] - other_configuration[i,:][2])**2 <= 0.49: # less than 0.7 angstrom
						dict_position[i]+=1
		list_conserve = list()
		for key in dict_position.keys():
			if dict_position[key] > 0.10*float(number_iteration):
				list_conserve.append(key)
				dict_count[list_vertices[key]]+=1
		
		for i in list_conserve:
			x_1, y_1, z_1 = configuration[i,:][0], configuration[i,:][1], configuration[i,:][2]
			for j in list_conserve:
				if j != i:
					x_2, y_2, z_2 = configuration[j,:][0], configuration[i,:][1], configuration[j,:][2]
					distance_local = (x_1 - x_2)**2 + (y_1 - y_2)**2 + (z_1 - z_2)**2
					distance_local = math.sqrt(distance_local)
					score = math.exp(-0.5*(distance_local- d[(list_vertices[i], list_vertices[j])])**2/0.5**2)
					random_prob = random.random()
					if score > random_prob:
						dict_corr[(list_vertices[i], list_vertices[j])].add(float('%.4f'%distance_local))
	return dict_corr, dict_count

def removekey(d, key):
	r = dict(d)
	del r[key]
	return r

number_iteration = 10
dict_model = iteration_procedure(d1, number_iteration)
dict_rmsd = dict()
dict_cluster = dict()
dict_state = dict()
for key in dict_model.keys():
	dict_state[key] = False
nb_cluster = 0

for key in dict_model.keys():
	if dict_state[key] == False: 
		nb_cluster+=1
		dict_cluster[nb_cluster] = set()
		dict_cluster[nb_cluster].add(key)
		configuration = dict_model[key]
		for key_other in dict_model.keys():
			other_configuration = dict_model[key_other]
			rmsd_local, other_configuration = kabsh(configuration, other_configuration)
			dict_rmsd[(key, key_other)] = rmsd_local
			print key, key_other, rmsd_local
			if rmsd_local <= 3.2:
				dict_state[key_other] = True
				dict_cluster[nb_cluster].add(key_other)
	print '__________'
	nb_rest = 0
	for key in dict_state.keys():
		if dict_state[key] == False:
			nb_rest +=1
	if nb_rest == len(dict_state.keys()):
		break

for key in dict_cluster.keys():
	if 1:#len(dict_cluster[key]) >= 20:
		print 'cluster: ', dict_cluster[key]
		print 'rmsd versus perfect cloud: '
		for element in dict_cluster[key]:
			configuration = dict_model[element]
			rmsd_, configuration = kabsh(configuration, D)
			print '%.4f'%rmsd_,
		print '\n'
		print '_______________'
			
	#rmsd_, configuration = kabsh(configuration, D)
	#dict_rmsd[(key, number_iteration + 1)] = rmsd_


def distance_dict_computation(configuration):
	d_ = dict()
	dict_point_ = dict()
	for i in range(len(list_vertices)):
		x_1, y_1, z_1 = configuration[i,:][0], configuration[i,:][1], configuration[i,:][2]
		dict_point_[list_vertices[i]] = Point(list_vertices[i], x_1, y_1, z_1)
	for key_1  in dict_point_.keys():
		for key_2 in dict_point_.keys():
			d_[(key_1, key_2)] = math.sqrt((dict_point_[key_1].x - dict_point_[key_2].x)**2+\
				(dict_point_[key_1].y - dict_point_[key_2].y)**2 +\
				(dict_point_[key_1].z - dict_point_[key_2].z)**2)
	return d_ 

def frobenius_norm(d_1, d_2):
	sum = 0.0
	for i in range(len(list_vertices)):
		for k in range(len(list_vertices)):
			sum+= (d_1[(list_vertices[i], list_vertices[k])] - d_2[(list_vertices[i], list_vertices[k])])**2
	return sum


















































#Influence=# tu as

















