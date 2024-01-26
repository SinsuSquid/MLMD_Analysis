import numpy as np
import matplotlib.pyplot as plt

plt.figure(figsize = (10,8))

data = np.genfromtxt('./dat/MSD.dat')

plt.plot(data[:,0], data[:,1],
         linewidth = 3, color = "#BE3455", label = "total")
plt.plot(data[:,0], data[:,2],
         linewidth = 3, color = "#6868AB", label = "Li$^+$")
plt.plot(data[:,0], data[:,3],
         linewidth = 3, color = "#939597", label = "Cl$^-$")
plt.plot(data[:,0], data[:,4],
         linewidth = 3, color = "#F5DF4D", label = "Al$^{3+}$")
plt.plot(data[:,0], data[:,0],
         linewidth = 3, color = 'k', linestyle = '--', label = "_nolegend_")

plt.title("Mean Squared Displacement")
plt.xlabel("Time (ps)")
plt.ylabel("MSD ($Å^2$)")

plt.xscale('log')
plt.yscale('log')

plt.xlim((0,1000))

plt.tight_layout()
plt.grid()
plt.legend()

plt.savefig('./figure/MSD.png', dpi = 300)
plt.show()
