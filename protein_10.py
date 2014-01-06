# increase 10% of all measurable distances
# compare created and original distance matrix to conclude.

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
  pass
