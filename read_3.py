# version 1.3
import re
import math
import numpy as np
import random
import itertools

def kabsh(P, Q):
  """ Kabsh algorithm to align two point clouds 
  return the rmsd between two structure and the aligned structure"""
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
  """ return the rmsd from two sets of vectors V and W """
  D = len(V[0])
  N = len(V)
  rmsd = 0.0
  for v, w in zip(V, W):
    rmsd+=sum((v[i]-w[i])**2 for i in range(D))
  return np.sqrt(rmsd/N)
  
def centroid(X):
  """return the centroid from a vectorset X """
  C = sum(X)/len(X)
  return C
  
class Point:
  def __init__(self, name = '', x=0.0, y=0.0, z=0.0):
    self.name = name
    self.x = x
    self.y = y
    self.z = z

dict_coord_complet = {}
file = open("/Users/James/protein_project/2KA0.pdb.1","r")
for line in file:
  line = line.rstrip('\n\r')
  line = re.split(" +", line)
  name = line[2] + '_' + line[5]
  if 'H' in name:
    coord_x = float(line[6])
    coord_y = float(line[7])
    coord_z = float(line[8])
    coord = Point(name, coord_x, coord_y, coord_z)
    dict_coord_complet[name] = point
  else:
    continue

vertices = set()
for item in dict_coord_complet.keys():
  vertices.add(item)

# generate the truncated and the complet distance matrix
full_matrix_file = open('/Users/James/protein_project/full_matrix_c','r')
distance_matrix_trunc = {}
distance_matrix_complet = {}
for line in full_matrix_file:
  line = line.rstrip('\n\r')
  line = re.split(' +', line)
  distance_matrix_complet[(line[0], line[1])] = float(line[2])
  if float(line[2]) <= 5.5:
    distance_matrix_trunc[(line[0], line[1])] = float(line[2])
  else:
    distance_matrix_trunc[(line[0], line[1])] = float('inf')

###################

def fw_algo(vertices, dict_distance):
  """run the FW algorithm to generate the temporary complet distance matrix from
  the truncated one """
  global vertices
  d = dict(dict_distance)
  for k in vertices:
    for i in vertices:
      for j in vertices:
        d[(i,j)] = min(d[(i,j)], d[(i,k)] + d[(k,j)])
  return d
  
###################
## construct the whole configuration from 4 points and a truncated distance matrix
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

####################################
D_complet = [ ] # configuration created from dict_coord_complet
list_vertices = list(vertices)
for key in list_vertices:
  x_D, y_D, z_D = dict_coord_complet[key].x, dict_coord_complet[key].y, dict_coord_complet[key].z
  coord_D = [x_D, y_D, z_D]
  D_complet.append(np.array(coord_D))
D_complet = np.array(D_complet)
Dc_complet = centroid(D_complet
D-=Dc_complet # turn the perfect cloud configuration around its centroid

#####################################

def iteration_procedure(d, number_iteration):
  global vertices, D_complet
  list_shuffle = list(vertices)
  dict_nb_model = {}
  dict_model = {} # gathering the contructed configurations
  
  nb_model = 0
  nb_try = 0
  while nb_try <= number_iteration:
    while 1:
      try: 
        random.shuffle(list_shuffle)
        list_4_points = list_shuffle[:4]
        dict_coord_constructed = constructing_4_points(list_4_points, distance_matrix_trunc)
        break
      except ValueError as err:
        pass
      except ZeroDivisionError as err:
        pass
    for i in range(4, len(list_shuffle)):
      next_point = list_shuffle[i]
      x_n, y_n, z_n = fifth_point(next_point, list_4_points, distance_matrix_trunc, dict_coord_constructed)
    
    V = [ ] # constructed configuration from the truncated distance matrix and the random 4 points
    for key in list(vertices):
      x_V, y_V, z_V = dict_coord_constructed[key].x, dict_coord_constructed[key].y, dict_coord_constructed[key].z
      coord_V = [x_V, y_V, z_V]
      V.append(np.array(coord_V))
    V = np.array(V)
    V_c = centroid(V)
    V-=V_c
    
    rmsd_kabsh_V, V = kabsh(V, D_complet)
    dict_model[nb_try] = V
    
    nb_try+=1
  return dict_model

def distance_matrix_generation(configuration):
  """ generating the distance matrix from a configuration """
  global vertices
  list_vertices = list(vertices)
  distance_matrix_gen = {}
  dict_constructed = dict()
  for i in range(len(list_vertices)):
    x_, y_, z_ = configuration[i,:][0], configuration[i,:][1], configuration[i,:][2]
    dict_constructed[list_vertice[i]] = Point(list_vertices[i], x_, y_, z_)
  for key_1 in dict_constructed.keys():
    for key_2 in dict_constructed.keys():
      distance_matrix_gen[(key_1, key_2)] = math.sqrt((dict_constructed[key_1].x - dict_constructed[key_2].x)**2+\
      (dict_constructed[key_1].y - dict_constructed[key_2].y)**2 + \
      (dict_constructed[key_1].z - dict_constructed[key_2].z)**2)
  return distance_matrix_gen
  
def harmonic_diff(d_1, d_2):
  """ returning the harmonic difference between two distance matrix
  especially against the initial distance matrix"""
  harmonic_diff = {}
  for key in d_1.keys():
    harmonic_diff[key] = (d_1[key] - d_2[key])**2
  return harmonic_diff
    
def frobenius_norm(harmonic_diff):
  """return the frobenius norm between two distance matrix """
  #harmonic_diff = harmonic_diff_matrix(d_1, d_2)
  sum = 0.0
  for key in harmonic_diff.keys():
    sum+= harmonic_diff[key]
  return math.sqrt(sum)/float(len(harmonic_diff.keys()))
  
##################  
def new_distance_matrix(dict_model):
  """finding consensus new distance matrix from different configurations """
  global vertices
  matrix_sum = {}
  for key_model in dict_model.keys():
  	configuration = dict_model[key_model]
  	distance_matrix_constructed = distance_matrix_generation(configuration)
  	harmonic_diff_matrix = harmonic_diff(distance_matrix_constructed, distance_matrix_trunc)
  	frobenius_norm_local = frobenius_norm(harmonic_diff_matrix)
  	for i in list_vertices:
  		for j in list_vertices:
  			matrix_sum[(list_vertices[i], list_vertices[j])] = \
  				distance_matrix_constructed[(list_vertices[i], list_vertices[j])]*\
  				
  pass
  
  
  
  
  
  
  
  
###################
# application

distance_matrix_trunc = fw_algo(vertices, distance_matrix_trunc) # complet the truncated matrix
number_iteration = 10
dict_model = iteration_procedure(distance_matrix_trunc, number_iteration)

  
  
  
  
  
  



















  
