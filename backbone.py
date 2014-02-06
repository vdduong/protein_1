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

