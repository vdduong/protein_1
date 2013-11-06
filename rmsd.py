"""
calculate rmsd between 2 XYZ files
by Jimmy Charnley Kromann
"""

import numpy
import sys
import re

def fit(P,Q):
  """
  varies the distance between P and Q, and optimizes rotation for each step until a minimum is found
  """
  
  step_size = P.max(0)
  threshold = step_size*1e-9
  rmsd_best = kabsch(P,Q)
  while True:
      for i in range(3):
        temp = numpy.zeros(3)
        temp[i] = step_size[i]
        rmsd_new = kabsch(P + temp, Q)
        if rmsd_new < rmsd_best:
            rmsd_best = rmsd_new
            P[:,i]+= step_size[i]
        else:
            rmsd_new = kabsch(P-temp, Q)
            if rmsd_new < rmsd_best:
              rmsd_best = rmsd_new
              P[:,i]-=step_size[i]
            else:
              step_size[i]/=2
        if (step_size < threshold).all():
          break
  return rmsd_best
  
  
def kabsh(P,Q):
  C = numpy.dot(numpy.transpose(P), Q)
  V, S, W = numpy.linalg.svd(C)
  d = (numpy.linalg.det(V)*numpy.linalg.det(W)) < 0.0
  if (d):
    S[-1] = - S[-1]
    V[:,-1] = - V[:,-1]
  
  U = numpy.dot(V,W)
  
  P = numpy.dot(P, U)
  return rmsd(P,Q)
  
def centroid(X):
  """ Calculate the centroid from a vectorset X """
  C = sum(X)/len(X)
  return C

def rmsd(V, W):
  """ calculate the rmsd from two sets of vectors V and W """
  D = len(V[0])
  N = len(V)
  rmsd = 0.0
  for v, w in zip(V, W):
    rmsd+=sum(([v[i]-w[i])**2.0 for i in range(D)])
  return numpy.sqrt(rmsd/N)
  
def get_coordinates(filename):
  f = open(filename, 'r')
  V = [ ]
  
  # skip 2 lines
  for _ in xrange(2):
    f.next()
  
  for line in f:
    numbers = re.findall(r'[-]?\d+\.\d+', line)
    numbers = [float(number) for number in numbers]
    V.append(numpy.array(numbers))
  f.close()
  V = numpy.array(V)
  return V

if __name__ == "__main__":
  args = sys.argv[1:]
  usage = ''' usage'''
  if len(args) < 2:
    print usage
    sys.exit(0)
    
  mol1 = args[0]
  mol2 = args[1]
  
  P = get_coordinates(mol1)
  Q = get_coordinates(mol2)
  
  print 'normal rmsd: ', rmsd(P,Q)
  Pc = centroid(P)
  Qc = centroid(Q)
  P-=Pc
  Q-=Qc
  print 'Kabsch rmsd: ', kabsch(P,Q)
  print 'fitted rmsd: ', fit(P,Q)
