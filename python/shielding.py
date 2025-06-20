#library to compute shielding functions by J.O. Wallace (1976) 
#see N-MISC-25-001
import numpy as np
import scipy.integrate as integrate

def sec(theta):
  return 1/np.cos(theta)

def csc(theta):
  return 1/np.sin(theta)

def cot(theta):
  return 1/np.tan(theta)

def slabFluxI(a,l1,l2,l3,mus=1.0):
  if(mus<=0): return 0

  integrand = lambda theta,phi: np.cos(theta)*(1-np.exp(-mus*l3*sec(phi)*sec(theta)))

  phiint = lambda phi: integrate.quad(integrand,0,np.math.atan(l2/((a+l3)*sec(phi))),args=(phi,))[0]
  thetint = integrate.quad(phiint,0,np.math.atan(l1/(a+l3)))[0]

  return thetint 

def slabFluxII(a,l1,l2,l3,mus=1.0):
  if(mus<=0): return 0

  integrand = lambda theta,phi: np.cos(theta)*(1-np.exp(-mus*(l1*cot(phi)-a)*sec(phi)*sec(theta)))

  phiint = lambda phi: integrate.quad(integrand,0,np.math.atan(l2/((a+l3)*sec(phi))),args=(phi,))[0]
  thetint = integrate.quad(phiint,np.math.atan(l1/(a+l3)),np.math.atan(l1/a))[0]

  return thetint 

def slabFluxIII(a,l1,l2,l3,mus=1.0):
  if(mus<=0): return 0

  integrand = lambda theta,phi: np.cos(theta)*(1-np.exp(-mus*(l2*cot(theta)-a)*sec(phi)*sec(theta)))

  phiint = lambda phi: integrate.quad(integrand,np.math.atan(l2/((a+l3)*sec(phi))),np.math.atan(l2/(a*sec(phi))),args=(phi,))[0]
  thetint = integrate.quad(phiint,0,np.math.atan(l1/(a+l3)))[0]

  return thetint 

def slabFluxIV(a,l1,l2,l3,mus=1.0):
  if(mus<=0): return 0

  integrand = lambda theta,phi: np.cos(theta)*(1-np.exp(-mus*(l1*cot(phi)-a)*sec(phi)*sec(theta)))

  phiint = lambda phi: integrate.quad(integrand,np.math.atan(l2/((a+l3)*sec(phi))),np.math.atan(l2/(l1*csc(phi))),args=(phi,))[0]
  thetint = integrate.quad(phiint,np.math.atan(l1/(a+l3)),np.math.atan(l1/a))[0]

  return thetint 

def slabFlux(a,l1,l2,l3,mus=1.0):
  return (1/(4*np.pi*mus))*(slabFluxI(a,l1,l2,l3,mus=mus)+slabFluxII(a,l1,l2,l3,mus=mus)+slabFluxIII(a,l1,l2,l3,mus=mus)+slabFluxIV(a,l1,l2,l3,mus=mus))
