import numpy as np
import matplotlib.pyplot as plt

vanHove = np.genfromtxt('./vanHove_distinct.out')
intermediate = np.genfromtxt('./intermediate_distinct.out')

plt.figure(figsize = (18,8))

plt.subplot(1,2,1)
cmap = plt.get_cmap('viridis', int((vanHove.shape[0] - 1) / 3))

x = vanHove[0,1:]
vanHove = vanHove[1:,1:]
for i in range(int(vanHove.shape[0] / 3)):
    plt.plot(x, vanHove[i,:], color = cmap(i), alpha = 0.2)

plt.xlabel("r (Å)")
plt.ylabel("$4 \pi r^2 G_d(r,t)$")
plt.xlim((0,10))

plt.grid()

plt.subplot(1,2,2)
for k in range(intermediate.shape[1] - 1):
    plt.plot(intermediate[1:,0], intermediate[1:,k+1],
             label = f'k = $2 \pi $ / {intermediate[0,k+1]:.2} ($Å^{-1}$)',
             markersize = 2, marker = 'o')

plt.xscale('log')
plt.ylim((-0.25, 1.25))
plt.xlim((1,1000))
plt.xlabel("Time (ps)")
plt.ylabel("$F_d(k,t)$")
plt.grid()
plt.legend()

plt.tight_layout()
plt.show()

