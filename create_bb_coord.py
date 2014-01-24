# compute backbone coordinates for HA, CA, HB, CB, HN, N
# for the protein 2KA0.pdb
import re
import math

file_out = open("/Users/James/protein_project/2KA0.bb","w")
file = open("/Users/James/protein_project/2KA0.pdb","r")
nb_line_read = 0
set_atom_authorized = ['N', 'CA', 'CB', 'HA', 'H', 'HB', 'HB2', 'HB3']
while nb_line_read <= 383:
  file.readline()
  nb_line_read+=1
for line in file:
  if 'ENDMDL' in line:
    break
  else:
    line = line.rstrip('\n\r')
    line = re.split(' +', line)
    if line[2] in set_atom_authorized:
      file_out.write('%s %s %s %s %s\n'%(line[2], line[5], line[6], line[7], line[8])) 
file_out.close()
file.close()
