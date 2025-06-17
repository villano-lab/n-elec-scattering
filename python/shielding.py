#library to compute shielding functions by J.O. Wallace (1976) 
#see N-MISC-25-001
import numpy as np

def sec(theta):
  return 1/np.cos(theta)

def slabFluxI(l1,l2,l3,mus):
  if(mus<=0): return 0

  integrand = lambda r,phi,theta: np.cos(theta)*(1-np.exp(-mus*l3*sec(phi)*sec(theta)))

  return integrand(5,0,0)
