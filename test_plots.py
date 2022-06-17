# %% IMPORTST 

import numpy as np
import matplotlib.pyplot as plt

# %% PARAMETERS AND FUNCTIONS

Up = 0.22
Ip = 0.5
omega = 0.057
N = 2
CEP = np.pi/2
ellip = 0

def envelope(t):
    return np.sin(omega * t / (2 * N))**2

def e_field(t):
    return 2 * np.sqrt(Up) * envelope(t) * omega * np.array([np.sin(omega * t + CEP) * np.cos(ellip / 2), -np.cos(omega * t + CEP) * np.sin(ellip / 2)])

def a_field(t):
    return 2 * np.sqrt(Up) * envelope(t) * np.array([np.cos(omega * t + CEP) * np.cos(ellip / 2), np.sin(omega * t + CEP) * np.sin(ellip / 2)])

