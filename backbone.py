# backbone

sd_helix = 4.3 # ??
sd_sheet = 2.7 # ??
sd_loop = 3.5 # ??
tol_distance = 0.05 # ??

def prob_helix(name_distance):
  return math.exp(-0.5*(dict_distance[name_distance]-sd_helix)**2/tol_distance**2)

def prob_sheet(name_distance):
  return math.exp(-0.5*(dict_distance[name_distance]-sd_sheet)**2/tol_distance**2)

def prob_loop(name_distance):
  return math.exp(-0.5*(dict_distance[name_distance]-sd_loop)**2/tol_distance**2)

# given the distance matrix between HN (correct or not), return the most probable sequential chain of HN

