#This is a small library that contains functions for the EDW theoretical backgrounds

import numpy as np
import pandas as pd
f3_cosmo = pd.read_csv('EDW/Fig3_cosmogenics.txt',comment="#",header=None)

#analytical forms of EDW III
def tritium(E):
    p0 = 1.406e-8
    p1 = 18.6
    p2 = 511
    return p0*(p1-E)**2*(p2+E)*(E**2+2*p2*E)**(1/2)

def surfbet(E):
    p0 = 1.34
    p1 = -0.058
    p2 = 0.2
    p3 = 40
    p4 = 11.4
    return p0*np.exp(p1*E)+p2*np.exp(-(E-p3)**2/(2*p4**2))

def lead(E):
    p0 = 0.037
    p1 = 0.15
    p2 = 95
    p3 = 5.7
    return p0+p1*np.exp(-(E-p2)**2/(2*p3**2))

def heat(E):
    p0 = 38.2725
    p1 = 0.293
    p2 = 1.4775
    p3 = 0.0812
    return p0*np.exp(-p1*E)+p2*np.exp(-p3*E)

def neutrons(E):
    p0 = 4.827e-4
    p1 = 0.3906
    p2 = 2.986e-4
    p3 = 0.05549
    return p0*np.exp(-p1*E)+p2*np.exp(-p3*E)

def cosmogenics(x):
    for i,y in enumerate(f3_cosmo[0]):
        if abs(x-y) < x/70: # 1% tolerance
            return f3_cosmo[1][i]
    return 0 #if nothing was returned

compton = 1e-1

def f3_total(x):
    return tritium(x)+surfbet(x)+lead(x)+heat(x)+neutrons(x)+[cosmogenics(val) for val in x]+compton

#print([cosmogenics(val) for val in f3_beta[0]]) #debug