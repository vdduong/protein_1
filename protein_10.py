# increase 10% of all measurable distances
# compare created and original distance matrix to conclude.
import numpy as np

def increase_10(d):
  for key_ in d.keys();
    if d[key_] < 5.5:
      d[key_] = d[key_]*1.1
    else:
      pass
  d = triangle_inequality(d)
  return d
  
def comparing_distance_matrix(d_new,d):
  sum = 0.0
  for key_ in d.keys():
    sum += d_new[key_]**2.0 - d[key_]**2.0
  sum = sum/2.0
  return sum
  
def update_distance_matrix():
  """ updating distance matrix from clusters of configurations and their respective difference against initial \
      distance matrix """
  global vertices
  # score based on distance matrix difference ???
  score = math.sqrt(comparing_distance_matrix(d_new, d))/len(vertices)
  
  # non linear scaled inter-atomic distance variance matrix
  sigma = np.zeros((len(vertices), len(vertices)))
  count_sigma = 0
  for configuration in dict_configuration:
    sigma[i][j]+= (d_new[i][j] - d[i][j])**2
    count_sigma+=1
  sigma = sigma/count_sigma
  
  # creating the matrix D
  D = np.zeros((len(vertices), len(vertices)))
  D[i][j] = norm((h[i] - h[j])) # norm to be defined
  # scaling non linear 
  h[i][j] = math.log(1 + sigma[i][j]**2/k_off)
  
  Z[diagonal] = 1 - 1/len(vertices)
  Z[non_diagonal] = - 1/len(vertices)
  
  G = -0.5*Z*D*Z.T
  
  
  break
