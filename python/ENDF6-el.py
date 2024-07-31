import ENDF6
import pandas as pds
import numpy as np

def fetch_elastic(filename='xn_data/si28_el.txt'):
  el = pds.read_csv(filename, skiprows=11,skipfooter=2, \
          names=['neutE', 'xn'],sep='\s+',engine='python')

  neute = np.asarray(el["neutE"],dtype=float)
  xn = np.asarray(el["xn"],dtype=float)

  return np.
