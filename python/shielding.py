#library to compute shielding functions by J.O. Wallace (1976) 
#see N-MISC-25-001
import pathlib
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parent
DATA_DIR = MODULE_DIR.parent / "data_files" 

import numpy as np
import ENDF6el as endfel
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

  phiint = lambda phi: integrate.quad(integrand,0,np.atan(l2/((a+l3)*sec(phi))),args=(phi,))[0]
  thetint = integrate.quad(phiint,0,np.atan(l1/(a+l3)))[0]

  return (1/(4*np.pi*mus))*thetint 

def slabFluxII(a,l1,l2,l3,mus=1.0):
  if(mus<=0): return 0

  integrand = lambda theta,phi: np.cos(theta)*(1-np.exp(-mus*(l1*cot(phi)-a)*sec(phi)*sec(theta)))

  phiint = lambda phi: integrate.quad(integrand,0,np.atan(l2/((a+l3)*sec(phi))),args=(phi,))[0]
  thetint = integrate.quad(phiint,np.atan(l1/(a+l3)),np.atan(l1/a))[0]

  return (1/(4*np.pi*mus))*thetint 

def slabFluxIII(a,l1,l2,l3,mus=1.0):
  if(mus<=0): return 0

  integrand = lambda theta,phi: np.cos(theta)*(1-np.exp(-mus*(l2*cot(theta)-a)*sec(phi)*sec(theta)))

  phiint = lambda phi: integrate.quad(integrand,np.atan(l2/((a+l3)*sec(phi))),np.atan(l2/(a*sec(phi))),args=(phi,))[0]
  thetint = integrate.quad(phiint,0,np.atan(l1/(a+l3)))[0]

  return (1/(4*np.pi*mus))*thetint 

def slabFluxIV(a,l1,l2,l3,mus=1.0):
  if(mus<=0): return 0

  integrand = lambda theta,phi: np.cos(theta)*(1-np.exp(-mus*(l1*cot(phi)-a)*sec(phi)*sec(theta)))

  phiint = lambda phi: integrate.quad(integrand,np.atan(l2/((a+l3)*sec(phi))),np.atan(l2/(l1*csc(phi))),args=(phi,))[0]
  thetint = integrate.quad(phiint,np.atan(l1/(a+l3)),np.atan(l1/a))[0]

  return (1/(4*np.pi*mus))*thetint 

def slabFlux(a,l1,l2,l3,mus=1.0):
  return (slabFluxI(a,l1,l2,l3,mus=mus)+slabFluxII(a,l1,l2,l3,mus=mus)+slabFluxIII(a,l1,l2,l3,mus=mus)+slabFluxIV(a,l1,l2,l3,mus=mus))

# get energies log scale
'''# Load total elastic cross sections for all 3 Silicon isotopes
f_28 = endfel.fetch_elastic(filename='../data_files/xn_data/si28_el.txt')
f_29 = endfel.fetch_elastic(filename='../data_files/xn_data/si29_el.txt')
f_30 = endfel.fetch_elastic(filename='../data_files/xn_data/si30_el.txt')

# Load total elastic cross sections for all 5 Germanium isotopes
f_70 = endfel.fetch_elastic(filename='../data_files/xn_data/ge70_el.txt')
f_72 = endfel.fetch_elastic(filename='../data_files/xn_data/ge72_el.txt')
f_73 = endfel.fetch_elastic(filename='../data_files/xn_data/ge73_el.txt')
f_74 = endfel.fetch_elastic(filename='../data_files/xn_data/ge74_el.txt')
f_76 = endfel.fetch_elastic(filename='../data_files/xn_data/ge76_el.txt')

ax1.plot(En, f_28(En), label="28-Si")
ax1.plot(En, f_29(En), label="29-Si", linestyle="--")
ax1.plot(En, f_30(En), label="30-Si", linestyle=":")
'''
    
def Edep28(a, l1, l2, l3, energies):
    f28=endfel.fetch_elastic(filename='../data_files/xn_data/si28_el.txt') #Load cross section from ENDF file for Si28 
    N=4.01501E23 # number density of oxygen
    flux=[]
    for i,E in enumerate(energies):
        sig=f28(E) #total elastic cross section
        big_sig=N*sig #macroscopic cross section
        flux.append(slabFlux(a, l1, l2, l3, big_sig))
    return flux
        
def Edep29(a, l1, l2, l3, energies):
    f29 = endfel.fetch_elastic(filename='../data_files/xn_data/si29_el.txt') #Load cross section from ENDF file for Si29
    N=4.01501E23 
    flux=[]
    for i,E in enumerate(energies):
        sig=f29(E) 
        big_sig=N*sig 
        flux.append(slabFlux(a, l1, l2, l3, big_sig))
    return flux



#Energy dependent flux functions for all five Germanium isotopes

def Edep70(a, l1, l2, l3, energies):
    f70 = endfel.fetch_elastic(filename='../data_files/xn_data/ge70_el.txt') #load cross section for Ge70
    N=4.01501E23 
    flux=[]
    for i,E in enumerate(energies):
        sig=f70(E) 
        big_sig=N*sig
        flux.append(slabFlux(a, l1, l2, l3, big_sig))
    return flux    


def Edep72(a, l1, l2, l3, energies):
    f72 = endfel.fetch_elastic(filename='../data_files/xn_data/ge72_el.txt') #Load cross section for Ge72
    N=4.01501E23 
    flux=[]
    for i,E in enumerate(energies):
        sig=f72(E)
        big_sig=N*sig 
        flux.append(slabFlux(a, l1, l2, l3, big_sig))
    return flux  

def Edep73(a, l1, l2, l3, energies):
    f73 = endfel.fetch_elastic(filename='../data_files/xn_data/ge73_el.txt') #Load cross section for Ge73
    N=4.01501E23 
    flux=[]
    for i,E in enumerate(energies):
        sig=f73(E)
        big_sig=N*sig 
        flux.append(slabFlux(a, l1, l2, l3, big_sig))
    return flux  

def Edep74(a, l1, l2, l3, energies):
    f74 = endfel.fetch_elastic(filename='../data_files/xn_data/ge74_el.txt')
    N=4.01501E23 
    flux=[]
    for i,E in enumerate(energies):
        sig=f74(E)
        big_sig=N*sig 
        flux.append(slabFlux(a, l1, l2, l3, big_sig))
    return flux  

def Edep76(a, l1, l2, l3, energies):
    f76 = endfel.fetch_elastic(filename='../data_files/xn_data/ge76_el.txt')
    N=4.01501E23 
    flux=[]
    for i,E in enumerate(energies):
        sig=f76(E)
        big_sig=N*sig 
        flux.append(slabFlux(a, l1, l2, l3, big_sig))
    return flux  

    
    #get total xn
    #run previous function (slab flux)
    #populate flux vector
    #return fluxvector
