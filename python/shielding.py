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


    
def Edep28(a, l1, l2, l3, energies):
    f28=endfel.fetch_elastic(filename=str(DATA_DIR/'xn_data'/'si28_el.txt')) #Load cross section from ENDF file for Si28 
    N=1.32941e23*0.922 # number density of Si in shotcrete * relative abundance of 28Si
    flux=[]
    for i,E in enumerate(energies):
        sig=f28(E) #total elastic cross section
        big_sig=N*sig*1e-24 #macroscopic cross section
        flux.append(slabFlux(a, l1, l2, l3, big_sig))
    return flux
        
def Edep29(a, l1, l2, l3, energies):
    f29 = endfel.fetch_elastic(filename=str(DATA_DIR/'xn_data'/'si29_el.txt')) #Load cross section from ENDF file for Si29
    N=1.32941e23*0.0466
    flux=[]
    for i,E in enumerate(energies):
        sig=f29(E) 
        big_sig=N*sig*1e-24 #number density * microscopic xn * barns to cm^2 conversion
        flux.append(slabFlux(a, l1, l2, l3, big_sig))
    return flux

def Edep30(a, l1, l2, l3, energies):
    f30 = endfel.fetch_elastic(filename=str(DATA_DIR/'xn_data'/'si30_el.txt')) #Load cross section from ENDF file for Si29
    N=1.32941e23*0.0307 
    flux=[]
    for i,E in enumerate(energies):
        sig=f30(E) 
        big_sig=N*sig*1e-24 #number density * microscopic xn * barns to cm^2 conversion
        flux.append(slabFlux(a, l1, l2, l3, big_sig))
    return flux

    

def Edep16(a, l1, l2, l3, energies):
    f16 = endfel.fetch_elastic(filename=str(DATA_DIR/'xn_data'/'o16_el.txt')) #Load cross section from ENDF file for O16
    N=4.01501E23 * 0.9975 #number density of O in shotcrete * relative abundance of 16O
    flux=[]
    for i,E in enumerate(energies):
        sig=f16(E) 
        big_sig=N*sig*1e-24  #number density of 16O * microscopic xn * barns to cm^2 conversion
        flux.append(slabFlux(a, l1, l2, l3, big_sig))
    return flux

def Edep17(a, l1, l2, l3, energies):
    f17 = endfel.fetch_elastic(filename=str(DATA_DIR/'xn_data'/'o17_el.txt')) #Load cross section from ENDF file for O17
    N=4.01501E23 *0.00038
    flux=[]
    for i,E in enumerate(energies):
        sig=f17(E) 
        big_sig=N*sig*1e-24
        flux.append(slabFlux(a, l1, l2, l3, big_sig))
    return flux

def Edep18(a, l1, l2, l3, energies):
    f18 = endfel.fetch_elastic(filename=str(DATA_DIR/'xn_data'/'o18_el.txt')) #Load cross section from ENDF file for O18
    N=4.01501E23 * 0.002
    flux=[]
    for i,E in enumerate(energies):
        sig=f18(E) 
        big_sig=N*sig*1e-24
        flux.append(slabFlux(a, l1, l2, l3, big_sig))
    return flux
    
def EdepTot(a, l1, l2, l3, energies):
    f28 = endfel.fetch_elastic(filename=str(DATA_DIR/'xn_data'/'si28_el.txt'))
    f29 = endfel.fetch_elastic(filename=str(DATA_DIR/'xn_data'/'si29_el.txt'))
    f30 = endfel.fetch_elastic(filename=str(DATA_DIR/'xn_data'/'si30_el.txt'))
    f16 = endfel.fetch_elastic(filename=str(DATA_DIR/'xn_data'/'o16_el.txt'))
    f17 = endfel.fetch_elastic(filename=str(DATA_DIR/'xn_data'/'o17_el.txt'))
    f18 = endfel.fetch_elastic(filename=str(DATA_DIR/'xn_data'/'o18_el.txt'))
    




    #return fluxvector
