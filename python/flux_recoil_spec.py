import numpy as np
import pandas as pd
import scipy.constants as co
import scipy.stats as ss
import itertools
import pickle
from scipy import signal
import ENDF6el as endfel
import masses as ms

# extrapolate line from lower-energy fast neutrons
E_thresh = 2e-2 # upper bound of linear region
E_therm = 0.15e-6 # near boundary of where thermal distribution has peak

def integrate_df(df):
    # (left-sided rectangular integral)
    dE = -df['E'].diff(periods = -1)
    dE.iat[-1] = dE.iat[-2]
    A = df['spec']*dE
    return A.sum()

def fast_extrapolation_line(E,fitted_line):
    # return height of fitted loglog line, if energy is larger than thermal threshold
    return np.exp(fitted_line.intercept + fitted_line.slope*np.log(E))*(E > E_therm)

def SNOLAB_flux_10keV(Emax=1):

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
  etotspec = np.zeros((np.size(etot),),dtype=np.float64)
  print(np.size(etot),np.size(etotspec))

  etotspec[(etot<=ffle[-1])] = fflespec_smooth
  etotspec[(etot>ffle[-1])] = ffhespec_smooth

  
  fast_lin_df = ffle[ffle < E_thresh]
  fast_lin_df_spec = fflespec_smooth[ffle< E_thresh]
  
  fitted_line = ss.linregress(np.log(fast_lin_df), np.log(fast_lin_df_spec))
  print(fitted_line)
  
  EE = np.geomspace(1000e-6, 2e-2, 10_000)

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

  #print(np.shape(En))
  dsig=np.zeros(np.shape(En))
  for i,E in enumerate(En):
    E*=1e6
    dsder = endfel.fetch_der_xn(En=E,Z=14,A=28,pts=1000,eps=1e-5)
    val = dsder(Er)
    if val>0:
      dsig[i] = val 
    #print(E,dsig[i])


  d = {'E':En,'spec':F*dsig}
  df = pd.DataFrame(data=d)
  #data=pd.DataFrame(np.array([En, F*dsig]), columns=['E', 'spec'])
  integral = integrate_df(df)

  return integral
