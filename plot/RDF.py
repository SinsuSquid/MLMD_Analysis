import numpy as np
import matplotlib.pyplot as plt

plt.figure(figsize = (10,8))

data = np.genfromtxt('./dat/RDF.dat')

plt.plot(data[:,0], data[:,1],
         linewidth = 5, color = "#BE3455", label = "Total")
plt.plot(data[:,0], data[:,2],
         linewidth = 2, color = "#6868AB", label = 'Li$^+$-Li$^+$',
         alpha = 0.5)
plt.plot(data[:,0], data[:,3],
         linewidth = 2, color = "#939597", label = 'Li$^+$-Cl$^-$',
         alpha = 0.5)
plt.plot(data[:,0], data[:,4],
         linewidth = 2, color = "#F5DF4D", label = 'Li$^+$-Al$^{3+}$',
         alpha = 0.5)
plt.plot(data[:,0], data[:,5],
         linewidth = 2, color = "#34568B", label = 'Cl$^-$-Cl$^-$',
         alpha = 0.5)
plt.plot(data[:,0], data[:,6],
         linewidth = 2, color = "#34568B", label = 'Cl$^-$-Al$^{3+}$',
         alpha = 0.5)
plt.plot(data[:,0], data[:,7],
         linewidth = 2, color = "#6B5B95", label = 'Al$^{3+}$-Al$^{3+}$',
         alpha = 0.5)

plt.title("Radial Distribution Function")

plt.xlim((0,8))

plt.xlabel("r (Å)")
plt.ylabel("RDF")

plt.tight_layout()
plt.grid()
plt.legend()

plt.savefig("./figure/RDF.png", dpi = 300)
plt.show()
