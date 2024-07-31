import ENDF6
import pandas as pds
import scipy
import numpy as np

def fetch_elastic(filename='xn_data/si28_el.txt'):
  el = pds.read_csv(filename, skiprows=11,skipfooter=2, \
          names=['neutE', 'xn'],sep='\s+',engine='python')

  neute = np.asarray(el["neutE"],dtype=float)
  xn = np.asarray(el["xn"],dtype=float)

  #make sure we are strictly increasing
  d = diff(neute)
  d=np.append(d,0)
  neute = neute[d>0]
  xn = xn[d>0]

  #get an interpolant function
  f=scipy.interpolate.UnivariateSpline(
        neute,
        xn,
        k=3,
        s=0,
        check_finite=True)

  return f 
