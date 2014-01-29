# obtaining C_alpha coordinates from proton clouds
# 
# algorithm: given the proton clouds coordinates
# calculate the values of internal distances and chirality
# approximate the coordinates of N, CO, C beta
#
# bond length CH = 1.09 angstrom
# bond length NH = 1.01
# bond length CC sp3 = 1.54
# bond length CN = 1.47

# the distance d_NN function of the backbone phi/psi angles
# d_NN = (13.6 + 2.69*cos(Phi) - 7.26*cos(Psi) - 1.26*cos(Phi)*cos(Psi) - 3.73*sin(Phi)*sin(Psi))*0.5
# in theory, d_NN between two adjacent amide protons can vary between the van der Waals limit of 2.0 and 
# a maximum distance of 4.7 angstrom 

import re
import random
import math

class Point:
  def __init__(self, name = '', x=0.0, y=0.0, z=0.0):
    self.name = name
    self.x = x
    self.y = y
    self.z = z


vertices = set()
vertices_HN = set()

dict_coord_complet = {}
file = open("/Users/James/protein_project/2KA0.pdb.1","r")
for line in file:
  line = line.rstrip('\n\r')
  line = re.split(" +", line)
  name = line[2] + '_' + line[5]
  if 'H' in name and 'O' not in name:
    coord_x = float(line[6])
    coord_y = float(line[7])
    coord_z = float(line[8])
    if line[2] == 'H':
  		vertices_HN.add(name)
    coord = Point(name, coord_x, coord_y, coord_z)
    dict_coord_complet[name] = coord
    
  else:
    continue


for item in dict_coord_complet.keys():
  if 'O' not in item:
	vertices.add(item)

print vertices_HN	

# generate the truncated and the complet distance matrix
full_matrix_file = open('/Users/James/protein_project/full_matrix_c','r')
distance_matrix_trunc = {}
distance_matrix_complet = {}
count_distance = 0
for line in full_matrix_file:
  line = line.rstrip('\n\r')
  line = re.split(' +', line)
  distance_matrix_complet[(line[0], line[1])] = float(line[2])
  if float(line[2]) <= 5.5:
  	count_distance +=1
  	#rand = random.random()
  	if count_distance%2 == 0:
  	  distance_matrix_trunc[(line[0], line[1])] = float(line[2])*1.0
  	else:
  		distance_matrix_trunc[(line[0], line[1])] = float(line[2])*4.0
  elif line[0] == line[1]:
  		distance_matrix_trunc[(line[0], line[1])] = 0.0	
  else:
    distance_matrix_trunc[(line[0], line[1])] = float('inf')

def fw_algo(dict_distance):
  """run the FW algorithm to generate the temporary complet distance matrix from
  the truncated one """
  global vertices
  d = dict(dict_distance)
  for k in vertices:
    for i in vertices:
      for j in vertices:
        d[(i,j)] = min(d[(i,j)], d[(i,k)] + d[(k,j)])
  return d
  
distance_matrix_trunc = fw_algo(distance_matrix_trunc)
distance_matrix_HH = dict()
for item_1 in vertices_HN:
	for item_2 in vertices_HN:
		distance_matrix_HH[(item_1, item_2)] = 0.0
for key in distance_matrix_HH.keys():
	distance_matrix_HH[key] = distance_matrix_trunc[key] # update

for key in distance_matrix_HH.keys():
	if 2.0 <=distance_matrix_HH[key] <=4.7:
		print key[0], key[1], distance_matrix_HH[key]		
# distance_matrix_NN
list_vertices_HN = list(vertices_HN)
while list_vertices_HN[0]!='H_20':
	random.shuffle(list_vertices_HN) 


def build_sequence_HN(list_vertices_HN):
	global distance_matrix_HH
	def total_distances(list_):
		distance = 0.0
		for i in xrange(0, len(list_)-1):
			distance+= distance_matrix_HH[(list_[i], list_[i+1])]
		return distance
		
	def compute_swap_indices(index, nb_HN):
		index_previous = (index - 1 + nb_HN)%nb_HN
		index_next = (index + 1)%nb_HN
		return (index_previous, index_next)
	
	def distance_swap(list_, index_a, index_b):
		index_A = min(index_a, index_b)
		index_B = max(index_a, index_b)
		
		(index_A_previous, index_A_next) = compute_swap_indices(index_A, len(list_))
		(index_B_previous, index_B_next) = compute_swap_indices(index_B, len(list_))
		
		distances = [ ]
		distances.append(distance_matrix_HH[(list_[index_A_previous], list_[index_A])])
		distances.append(distance_matrix_HH[(list_[index_B], list_[index_B_next])])
		if index_A == index_B_previous:
			distances.append(distance_matrix_HH[(list_[index_A], list_[index_B])])
		else:
			distances.append(distance_matrix_HH[(list_[index_A], list_[index_A_next])])
			distances.append(distance_matrix_HH[(list_[index_B_previous], list_[index_B])])
		return sum(distances)
		
	def annealing(list_vertices_HN, temperature_begin=1.0e+300, temperature_end=0.1, cooling_factor=.999, nb_iterations=1):
		list_best = list_vertices_HN[:]
		distance_best = total_distances(list_best)
		
		distances_current = []
		distances_best = []
		ids_iteration = []
		
		try:
			for iteration in range(nb_iterations):
				temperature = temperature_begin
				list_current = list_best[:]
				distance_current = distance_best
				distance_new = distance_best
				list_new = list_best[:]
				
				step = 0
				while temperature > temperature_end:
					index = random.sample(range(len(list_new)-1),2)
					index[0] +=1
					index[1] +=1
					
					swap_before = distance_swap(list_new, index[0], index[1])
					list_new[index[0]], list_new[index[1]] = list_new[index[1]], list_new[index[0]]
					swap_after = distance_swap(list_new, index[0], index[1])
					
					# compute new distance
					distance_new = distance_new - swap_before + swap_after
					
					diff = distance_new - distance_current
					if diff < 0 or math.exp(-diff/temperature) > random.random():
						list_current = list_new[:]
						distance_current = distance_new
					else:
						distance_new = distance_current
						list_new = list_current[:]
					
					if distance_current < distance_best:
						list_best = list_current[:]
						distance_best = distance_current
					print list_best, distance_best
					
					if True:
						distances_current.append(distance_current)
						distances_best.append(distance_best)
					temperature = temperature*cooling_factor
					step+=1					
					
		except KeyboardInterrupt, e:
			print 'interrupt'
		return list_best, distance_best
	list_best, distance_best = annealing(list_vertices_HN, temperature_begin=1.0e+300, temperature_end=0.1, cooling_factor=.99, nb_iterations=1)
	return list_best, distance_best

#list_best, distance_best = build_sequence_HN(list_vertices_HN)
#print '_____'
#print list_best
#print distance_best
#def build_NN_backbone():
#	if 2.0 <= distance_matrix_NN[(list_vertices_N[i], list_vertices_N[j])] <= 4.7:





			
