import sys
import numpy as np
from scipy.stats import linregress

data = np.genfromtxt(sys.argv[1])

x = data[:,0].flatten()
li_msd = data[:,1].flatten()

regress = linregress(x[200:1000], li_msd[200:1000])

D = regress.slope / 6.0

print(f"D = {D} Å^2 / ps")
A2cm = 1.0E-8
ps2s = 1.0E-12

D = D * A2cm * A2cm
D = D / ps2s

print(f"D = {D} cm^2 / s")

k = 8.6173E-5 # ev / K
T = 300 # K
n = 432
n /= 42.123986964 * A2cm
n /= 39.1151307524 * A2cm
n /= 39.1151307524 * A2cm # / cm^3
e = 1.60217663E-19 # C
z = 1

sigma = (n * e * z * z) / (k * T) * D # A / V cm^-1 (S cm^-1) 

print(f"sigma = {sigma:e} S / cm")
