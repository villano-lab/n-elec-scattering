import numpy as np
import pandas as pd
import math
import random

# Define Isotope class
class Isotope:
    def __init__(self, name, half_life, decay_modes, count=0, is_stable=False):

        self.name = name
        self.half_life = half_life
        self.decay_modes = decay_modes
        self.count = count               
        self.is_stable = is_stable 
    
    # Default "count" is zero, set to different count for OG parent isotopes
    # "decay_modes" is list of tuples with decay mode (alpha or beta), daughter, and branching ratio
    
    def decay(self):
        if self.is_stable or not self.decay_modes:
            return ("stable", self.name)

        r = random.random()
        cumulative = 0.0
        for mode, daughter, prob in self.decay_modes:
            cumulative += prob
            if r < cumulative:
                return (mode, daughter)
        # Fallback in case of rounding errors
        return self.decay_modes[-1]

    def time_until_decay(self):
        return random.expovariate(math.log(2) / self.half_life)

    
## Thorium decay chain isotope dictionary
Th_decay_data = {
    "Ra_228": {
        "half_life": 5.7 * 365 * 24 * 3600,
        "decay_modes": [("beta", "Ac_228", 1.0)]
    },
    "Ac_228": {
        "half_life": 6.1 * 60,
        "decay_modes": [("beta", "Th_228", 1.0)]
    },
    "Th_228": {
        "half_life": 1.9 * 365 * 24 * 3600,
        "decay_modes": [("alpha", "Ra_224", 1.0)]
    },
    "Ra_224": {
        "half_life": 3.6 * 24 * 3600,
        "decay_modes": [("alpha", "Rn_220", 1.0)]
    },
    "Rn_220": {
        "half_life": 55,
        "decay_modes": [("alpha", "Po_216", 1.0)]
    },
    "Po_216":{
        "half_life":  0.14,
        "decay_modes": [("alpha", "Pb_212", 1.0)]
    },
    "Pb_212":{
        "half_life":  10.6 * 3600,
        "decay_modes": [("beta", "Bi_212", 1.0)]
    },
    "Bi_212":{
        "half_life":  61 * 60,
        "decay_modes": [("beta", "Po_212", 0.6406),
                        ("alpha", "Tl_208", 0.3594)]
    },
    "Po_212":{
        "half_life":  3e-7,
        "decay_modes": [("alpha","Pb_208", 1.0)]
    },
    "Tl_208":{
        "half_life": 3.1 * 60,
        "decay_modes": [("beta", "Pb_208", 1.0)]}
}

#Define OG parent, identify count, and define first stable daughter (end of chain)
Th_232 = Isotope("Th_232", 1.41e10 * 365 * 24 * 3600, [("alpha", "Ra_228", 1.0)], 6e23)
Pb_208 = Isotope("Pb_208", float('inf'), [], is_stable=True)

# Build Isotope objects
th_chain_isotopes = {}
for name, data in Th_decay_data.items():
    th_chain_isotopes[name] = Isotope(name, data["half_life"], data["decay_modes"])

# List of thorium chain Isotope instances
th_chain_isotope_list = [Th_232] + list(th_chain_isotopes.values()) + [Pb_208]


## Uranium decay chain dictionary
U_decay_data = {
    "Th_231": {
        "half_life": 25.52 * 3600,
        "decay_modes": [("beta", "Pa_231", 1.0)]
    },
    "Pa_231": {
        "half_life": 32760 * 365 * 24 * 3600,
        "decay_modes": [("beta", "Ac_227", 1.0)]
    },
    "Ac_227": {
        "half_life": 21.772 * 365 * 24 * 3600,
        "decay_modes": [("alpha", "Fr_223", 0.0138),
                       ("beta", "Th_227", 0.9862)]
    },
    "Th_227": {
        "half_life": 18.68 *24*3600,
        "decay_modes": [("alpha", "Ra_223", 1)]
    },
    "Ra_223": {
        "half_life": 11.43 *24*3600,
        "decay_modes": [("alpha", "Rn_219", 1)]
    },
    "Fr_223": {
        "half_life": 22 * 60,
        "decay_modes": [("alpha", "Ra_223", 0.00006),
                       ("beta", "At_219", 0.99994)]
    },
    "At_219": {
        "half_life": 56,
        "decay_modes": [("alpha", "Bi_215", 0.936),
                       ("beta", "Rn_219", 0.064)]
    },
    "Rn_219": {
        "half_life": 3.93,
        "decay_modes": [("alpha", "Po_215", 1)]
    },
    "Bi_215": {
        "half_life": 7.6 * 60,
        "decay_modes": [("beta", "Po_215", 1)]
    },
    "Po_215": {
        "half_life": 1.781e-3,
        "decay_modes": [("alpha", "Pb_211", 0.9999977),
                       ("beta", "At_215", 0.0000023)]
    },
    "At_215": {
        "half_life": 1e-4,
        "decay_modes": [("alpha", "Bi_211", 1)]
    },
    "Pb_211":{
        "half_life": 36.1 * 60,
        "decay_modes": [("beta", "Bi_211", 1)]
    },
    "Bi_211":{
        "half_life": 2.14 * 60,
        "decay_modes": [("alpha", "Tl_207", 0.99724), 
                        ("beta","Po_211", 0.00276)]
    },
    "Po_211":{
        "half_life": 0.516,
        "decay_modes":[("alpha", "Pb_207", 1)]
    },
    "Tl_207":{
        "half_life": 4.77 * 60,
        "decay_modes":[("beta", "Pb_207", 1)]
    }
}

# Define OG parent outside dictionary, identify count, and first stable daughter (end of chain)
U_235 = Isotope("U_235", 7.04e8*365*24*3600, [("alpha", "Th_231", 1.0)], 6e23)
Pb_207 = Isotope("Pb_207", float('inf'), [], is_stable=True)

# Build Isotope objects for uranium decay chain
u_chain_isotopes = {}
for name, data in U_decay_data.items():
    u_chain_isotopes[name] = Isotope(name, data["half_life"], data["decay_modes"])

# List of Isotope instances
u_chain_isotope_list = [U_235] + list(u_chain_isotopes.values()) + [Pb_207]
    

    
