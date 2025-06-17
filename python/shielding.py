#library to compute shielding functions by J.O. Wallace (1976) 
#see N-MISC-25-001
import numpy as np
import scipy.integrate as integrate

def sec(theta):
  return 1/np.cos(theta)

def slabFluxI(a,l1,l2,l3,mus=1.0):
  if(mus<=0): return 0

  integrand = lambda theta,phi: np.cos(theta)*(1-np.exp(-mus*l3*sec(phi)*sec(theta)))

  phiint = lambda phi: integrate.quad(integrand,0,np.math.atan(l2/(a+l3)*sec(phi)),args=(phi,))[0]
  thetint = integrate.quad(phiint,0,np.math.atan(l1/(a+l3)))[0]

  return thetint 
