# coding: utf-8

# # APN Critical aggregation model
# # "APN"-model for the concentrations of monomers and aggregates near the critical aggregation concentration
# 
# More details on our web page
# http://www.usc.es/fotofqm/en/units/single-molecule-fluorescence/concentration-model-surfactants-near-cmc
# 
# A Model for Monomer and Micellar Concentrations in Surfactant Solutions. Application to Conductivity, NMR, Diffusion and Surface Tension data
# Wajih Al-Soufi, Lucas Piñeiro, Mercedes Novo, Journal of Colloid and Interface Science 2012, 370, 102–110 DOI: 10.1016/j.jcis.2011.12.037
#  
# 
# V2 - 2023/10/31

import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.special as sc

# ### Concentration Model
# #### APNS1(cS0,cmc,r)
# Concentration model, used by all other derived properties.
# Takes [S]0 and calculates the monomeric concentration [S1] as function of the cmc and the relative transition width r 

def APNS1(cS0,cmc,r):
    s0 = cS0/cmc;
    pi = math.pi;
    A=2/(1+np.sqrt(2/pi)*r*np.exp(-1/(2*r*r))+sc.erf(1/np.sqrt(2)/r));
    cS1=cmc*(1-(A/2)*(np.sqrt(2/pi)*r*np.exp(-(s0-1)**2/(2*r*r))+(s0-1)*(sc.erf((s0-1)/(np.sqrt(2)*r))-1)));    
    return cS1


# ### Direct Surfactant Properties:
# #### APNConductivity(cS0,cmc,r,a,b,c)
# Takes [S]0 and calculates the electric conductivity κ as function of the cmc, the relative transition width r, 
# the slopes a and b, and the solvent conductivity c = κs.

def APNConductivity(cS0,cmc,r,a,b,c):
    cS1 = APNS1(cS0, cmc, r)
    cSm = cS0 - cS1       
    return a * cS1 + b * cSm + c
