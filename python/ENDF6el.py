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
  d = np.diff(neute)
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

def fetch_elastic_angular(filename='xn_data/n-014_Si_028.endf'):

  #check if it's an ENDF file somehow?

  #get the section of the ENDF file that has ang dists
  f = open(filename)
  lines = f.readlines()
  sec = ENDF6.find_section(lines, MF=4, MT=2)

  #check if it has a legendre coeffs section?--I'm assuming it does right now..

  #get number of data points
  arr=np.str_.split(sec[3])
  num=int(arr[0])
  print(num)
  al = np.zeros((num,36)) #assume no more than 36 legendre coeffs
  en = np.zeros((num,))

  #get the coeffs into the data structure
  readl=True
  linecnt=0
  ecnt=0
  tote=num
  mpoles=0
  for ln in sec[4:-1]:
      #print(ln)
      #break away if you're done reading e points
      if ecnt==tote:
        break
      #read the number of multipoles
      if readl:
        arr=np.str_.split(ln)
        mpoles=int(arr[4])
        en[ecnt]=float(arr[1].replace("-", "e-").replace("+", "e+").lstrip("e"))
        #print(mpoles)
        readl=False
      else:
        #arr=np.str_.split(ln)
        a = ENDF6.read_line(ln)
        #print(a)
        i1=linecnt*6
        i2=i1+6
        #print(i1,i2)
        al[ecnt,i1:i2] = a
        #ecnt+=1
        linecnt+=1
        if linecnt==np.ceil(mpoles/6.0):
           linecnt=0
           ecnt+=1
           readl=True
          
      arr = np.char.split(ln)
      #print(np.shape(arr))
      if(np.shape(arr)==()): continue
      #print([np.fromstring(x) for x in arr])
 
  #close file
  f.close() 

  return (en,al)



