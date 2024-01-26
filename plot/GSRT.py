import numpy as np
import matplotlib.pyplot as plt

plt.figure(figsize = (10,8))
data = np.genfromtxt('./dat/GSRT.dat')

plt.plot(data[:,0], data[:,1],
         linewidth = 3, color = "#BE3455", label = "100 ps")
plt.plot(data[:,0], data[:,2],
         linewidth = 3, color = "#6868AB", label = "500 ps")
plt.plot(data[:,0], data[:,3],
         linewidth = 3, color = "#939597", label = "1000 ps")
plt.plot(data[:,0], data[:,4],
         linewidth = 3, color = "#F5DF4D", label = "1500 ps")
plt.plot(data[:,0], data[:,5],
         linewidth = 3, color = "#34568B", label = "2000 ps")

plt.xlim((0, 15))

plt.vlines(4.323, 0.000, 0.010, linestyle = '--', color = 'k')

plt.title("van Hove correlation of Li$^+$ ion")
plt.xlabel("r (Å)")
plt.ylabel("4πr$^2$ $G_s(r,t)$")

plt.tight_layout()
plt.legend()
plt.grid()

plt.savefig('./figure/GSRT.png', dpi = 300)
plt.show()
