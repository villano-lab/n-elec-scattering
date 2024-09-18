import numpy as np
import pandas as pd
import scipy.constants as co
import scipy.stats as ss
import itertools
import pickle
from scipy import signal

def integrate_df(df):
    # (left-sided rectangular integral)
    dE = -df['E'].diff(periods = -1)
    dE.iat[-1] = dE.iat[-2]
    A = df['spec']*dE
    return A.sum()

def SNOLAB_flux(Emax=1):

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

  return 
