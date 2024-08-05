import scipy.constants as co
import fortranformat as ff
import numpy as np

mnbar = 931.494045;
m_n = co.physical_constants['neutron mass energy equivalent in MeV'][0]
m_e = co.physical_constants['electron mass energy equivalent in MeV'][0]


def readFile(filename='isotope_data/mass_1.mas20.txt'):

  format = ff.FortranRecordReader('(a1,i3,i5,i5,i5,1x,a3,a4,1x,f14.6,f12.6,f13.5,1x,f10.5,1x,a2,f13.5,f11.5,1x,i3,1x,f13.6,f12.6)')
  f = open(filename, 'r')
  masses = []
  count = 0
  for line in f:
      count = count + 1
      if count > 36:
          masses.append(format.read(line.replace('*', ' ').replace('#', '.')))
  f.close()
  return masses


def getMass(Z0=0,A0=1):

  masses = readFile()
  Z=np.zeros((len(masses),))
  A=np.zeros((len(masses),))
  delta=np.zeros((len(masses),))
  edelta=np.zeros((len(masses),))
  
  for i,a in enumerate(masses):
    delta[i] = masses[i][7]
    edelta[i] = masses[i][7]
    Z[i] = masses[i][3]
    A[i] = masses[i][4]

  azdel = delta[(Z==Z0)&(A==A0)]

  return (A0)*mnbar + (azdel/1000) -Z0*m_e;
