import numpy as np
import pandas as pd
import scipy.constants as co
import scipy.stats as ss
import itertools
import pickle
from scipy import signal
import ENDF6el as endfel
import masses as ms
import scipy.constants as co
import periodictable as pt

#############
# Constants #
#############

#constants
NA = co.physical_constants['Avogadro constant'][0]
s2day = 1/(60*60*24) #seconds to days
#constants for calcs, first in SI units
gn = co.physical_constants['neutron gyromag. ratio'][0] #default is s^-1 T^-1; CGS is s^-1 Gauss^-1
mub = co.physical_constants['Bohr magneton'][0] #default is J T^-1
hbar = co.physical_constants['reduced Planck constant'][0] #default in J s

#convert to CGS
#see https://en.wikipedia.org/wiki/Centimetre%E2%80%93gram%E2%80%93second_system_of_units
m_n_CGS = co.physical_constants['neutron mass'][0]*1e3 #convert to grams
gn_CGS = gn/1e4
mub_CGS = mub*1e3
hbar_CGS = hbar*1e7

# extrapolate line from lower-energy fast neutrons
E_thresh = 2e-2 # upper bound of linear region
E_therm = 0.15e-6 # near boundary of where thermal distribution has peak

#############
# Functions #
#############

def integrate_df(df):
    # (left-sided rectangular integral)
    dE = -df['E'].diff(periods = -1)
    dE.iat[-1] = dE.iat[-2]
    A = df['spec']*dE
    return A.sum()

def fast_extrapolation_line(E,fitted_line):
    # return height of fitted loglog line, if energy is larger than thermal threshold
    return np.exp(fitted_line.intercept + fitted_line.slope*np.log(E))*(E > E_therm)

def Emax(En): #En in keV; returns maximum recoil energy for neutron energy
    return (4*ms.m_e*ms.m_n*En)/(ms.m_e+ms.m_n)**2

def Enmin(Er): #recoil energy in keV; returns minimum neutron energy to give that recoil energy
    return (Er*(ms.m_e+ms.m_n)**2)/(4*ms.m_e*ms.m_n)

def dsigdErNE(En,Er):
    if(Er<Emax(En)):
      return 8*np.pi*m_n_CGS**2*gn_CGS**2*mub_CGS**2/hbar_CGS**2/Emax(En)
    else: 
      return 0

def SNOLAB_flux(Enmin=1e-3):

  # read in fast neutron flux spectrum (from reading_n_spectra.ipynb)
  fast_flux_df = pd.read_pickle('../data_files/FDF.txt') # 'E' in MeV, 'spec' in neutrons cm^-2 sec^-1 MeV^-1

  #use numpy arrays
  ff = np.asarray(fast_flux_df['E']);
  ffspec = np.asarray(fast_flux_df['spec']);

  # calculate flux level of fast neutrons
  fast_flux = integrate_df(fast_flux_df)
  print('fast flux: {} n/m^2/day'.format((fast_flux*10000)*24*60*60))

  #smooth the data
  ffspec_smooth = signal.savgol_filter(ffspec, 501, 3) # window size 501, polynomial order 3

  cutoff=0.3

  ffhe = ff[ff>cutoff]
  ffhespec = ffspec[ff>cutoff]
  
  #smooth the data
  ffhespec_smooth = signal.savgol_filter(ffhespec, 2001, 3) # window size 1001, polynomial order 3
  
  ffle = ff[ff<=cutoff]
  fflespec = ffspec[ff<=cutoff]
  print(np.size(ffle))
  
  #smooth the data
  fflespec_smooth = signal.savgol_filter(fflespec, 75, 3) # window size 1001, polynomial order 3

  etot = np.concatenate((ffle,ffhe))
  etot = np.unique(etot)
  print('shape of etot: {}'.format(np.shape(etot)))
  etotspec = np.zeros((np.size(etot),),dtype=np.float64)
  print(np.size(etot),np.size(etotspec))

  etotspec[(etot<=ffle[-1])] = fflespec_smooth
  etotspec[(etot>ffle[-1])] = ffhespec_smooth

  
  fast_lin_df = ffle[ffle < E_thresh]
  fast_lin_df_spec = fflespec_smooth[ffle< E_thresh]
  
  fitted_line = ss.linregress(np.log(fast_lin_df), np.log(fast_lin_df_spec))
  print(fitted_line)
  
  EE = np.geomspace(Enmin, 2e-2, 10_000)

  #ax1.plot(ff, ffspec,label='simulated flux')
  #ax1.plot(etot, etotspec,color='orange',label="smoothed flux")
  #plt.plot(EE, fast_extrapolation_line(EE), color='orange', linestyle = 'dashed',label="extrapolated flux")

  #trim the extrapolated part
  upperE=np.min(etot)
  EE=EE[EE<upperE]
  EF=fast_extrapolation_line(EE,fitted_line)

  #cat them together
  E=np.append(EE,etot)
  F=np.append(EF,etotspec)

  print(np.max(EE),np.min(etot))

  return E,F,ff,ffspec

def dRdEr(Er,En,F,N=100,Z=14,A=28):


  #get min neutron energy
  mass = ms.getMass(Z,A)
  Enmin = Er/(4*mass*ms.m_n/(mass+ms.m_n)**2)
  #print(Enmin)

  #cut down the density of points
  idx=np.arange(0,len(En),1)
  cidx=idx%N==0
  En=En[cidx]
  F=F[cidx]

  #trim the flux energies for ones that can actually contribute
  cEn=En>=Enmin
  En=En[cEn]
  F=F[cEn]

  if(np.shape(En)[0]<2):
    return 0.0


  #get the filenames
  symbol=pt.elements[Z].symbol
  symbol_lower=symbol.lower()
  sigtotfile='../data_files/xn_data/{0:}{1:}_el.txt'.format(symbol_lower,A)
  endffile='../data_files/xn_data/n-{0:03d}_{1:}_{2:03d}.endf'.format(Z,symbol,A)
  print(sigtotfile,endffile)


  #print(np.shape(En))
  dsig=np.zeros(np.shape(En))
  for i,E in enumerate(En):
    E*=1e6
    dsder = endfel.fetch_der_xn(En=E,M=mass,pts=1000,eps=1e-5,sigtotfile=sigtotfile,endffile=endffile)
    val = dsder(Er)
    if val>0:
      dsig[i] = val 
    #print(E,dsig[i])


  dsig*=(endfel.barns2cm2*endfel.keV2MeV)
  d = {'E':En,'spec':F*dsig}
  df = pd.DataFrame(data=d)
  #data=pd.DataFrame(np.array([En, F*dsig]), columns=['E', 'spec'])
  integral = integrate_df(df)

  integral*=(1/s2day) 
  integral*=(1/(ms.getAMU(Z,A)*ms.amu2g*1e-3))  
  return integral

def dRdErfast(Er,En,F,N=100,Z=14,A=28):

  #vectorize Er
  if isinstance(Er, float):
        Er=[Er]
  Er = np.asarray(Er)
  #print(Er)

  #get min neutron energy
  mass = ms.getMass(Z,A)
  Enmin = Er/(4*mass*ms.m_n/(mass+ms.m_n)**2)
  #print(Enmin)

  #cut down the density of points
  idx=np.arange(0,len(En),1)
  cidx=idx%N==0
  En=En[cidx]
  F=F[cidx]

  if(np.shape(En)[0]<2):
    return 0.0

  #get the filenames
  #sigtotfile='../data_files/xn_data/si{}_el.txt'.format(A)
  #endffile='../data_files/xn_data/n-{0:03d}_Si_{1:03d}.endf'.format(Z,A)
  #print(sigtotfile,endffile)

  #get the filenames
  symbol=pt.elements[Z].symbol
  symbol_lower=symbol.lower()
  sigtotfile='../data_files/xn_data/{0:}{1:}_el.txt'.format(symbol_lower,A)
  endffile='../data_files/xn_data/n-{0:03d}_{1:}_{2:03d}.endf'.format(Z,symbol,A)
  print(sigtotfile,endffile)

  #make big ole matrix
  dsig=np.zeros((np.shape(En)[0],np.shape(Er)[0]))
  for i,E in enumerate(En):
    E*=1e6
    dsder = endfel.fetch_der_xn(En=E,M=mass,pts=1000,eps=1e-5,sigtotfile=sigtotfile,endffile=endffile)
    dsig[i,:] = dsder(Er)
    #print(dsig[i,:])


  #remove negatives 
  dsig[dsig<0]=0


  #integrate
  integral=np.zeros(np.shape(Er))
  for i,E in enumerate(Er): 

    #trim the flux energies for ones that can actually contribute
    enidx=np.arange(0,len(En),1)
    cEn=En>=Enmin[i]
    Entemp=En[cEn]
    Ftemp=F[cEn]
    enidx=enidx[cEn]
    xn=dsig[enidx,i]
    xn*=(endfel.barns2cm2*endfel.keV2MeV)
    if(np.shape(enidx)[0]<2):
      integral[i]=-999999999
    else:
      #print(np.shape(enidx))
      d = {'E':Entemp,'spec':Ftemp*xn}
      df = pd.DataFrame(data=d)
      #data=pd.DataFrame(np.array([En, F*dsig]), columns=['E', 'spec'])
      integral[i] = integrate_df(df)


  integral*=(1/s2day) 
  integral*=(1/(ms.getAMU(Z,A)*ms.amu2g*1e-3))  
  return integral,dsig

def dRdErNE(Er,En,F,N=1,Z=14,A=28,eta=1): #for neutron scattering off electrons eta is the electrons/atom

  #vectorize Er
  if isinstance(Er, float):
        Er=[Er]
  Er = np.asarray(Er)
  #print(Er)

  #get min neutron energy
  Enmin = Er/(4*ms.m_e*ms.m_n/(ms.m_e+ms.m_n)**2)
  #print(Enmin)

  #cut down the density of points
  idx=np.arange(0,len(En),1)
  cidx=idx%N==0
  En=En[cidx]
  F=F[cidx]

  if(np.shape(En)[0]<2):
    return 0.0

  #make big ole matrix
  dsig=np.zeros((np.shape(En)[0],np.shape(Er)[0]))
  for i,E in enumerate(En):
    E*=1e3 #E in keV
    dsder = lambda x: dsigdErNE(E,x)
    dsderv = np.vectorize(dsder)
    dsig[i,:] = dsderv(Er*1e3)
    #print(E,dsig[i,:])


  #remove negatives 
  dsig[dsig<0]=0


  #integrate
  integral=np.zeros(np.shape(Er))
  for i,E in enumerate(Er): 

    #trim the flux energies for ones that can actually contribute
    enidx=np.arange(0,len(En),1)
    cEn=En>=Enmin[i]
    Entemp=En[cEn]
    Ftemp=F[cEn]
    enidx=enidx[cEn]
    xn=dsig[enidx,i]
    #xn*=(endfel.barns2cm2*endfel.keV2MeV)
    #xn*=(endfel.barns2cm2)
    if(np.shape(enidx)[0]<2):
      integral[i]=-999999999
    else:
      #print(np.shape(enidx))
      d = {'E':Entemp,'spec':Ftemp*xn}
      df = pd.DataFrame(data=d)
      #data=pd.DataFrame(np.array([En, F*dsig]), columns=['E', 'spec'])
      integral[i] = integrate_df(df)


  integral*=(1/s2day) 
  integral*=(eta/(ms.getAMU(Z,A)*ms.amu2g*1e-3))  
  return integral,dsig
