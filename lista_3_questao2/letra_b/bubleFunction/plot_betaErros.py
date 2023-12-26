import numpy as np
import os
import matplotlib.pyplot as plt
import math

file = 'errorL2.csv'

x_projL2 = np.genfromtxt(file, delimiter=',', max_rows=1)
u_projL2 = np.genfromtxt(file, delimiter=',', skip_header=1, max_rows=1)

minError_idx = np.argmin(np.array(u_projL2))
beta_minError = x_projL2[minError_idx]
print(r'$\beta$'," = ", beta_minError)
print("MinError = ", u_projL2[minError_idx])


plt.figure()
plt.grid(True)

plt.plot(x_projL2, u_projL2, c="r", marker='o', linestyle='--', label=r'$\beta$', zorder=2)
    
plt.xlabel("x")
plt.ylabel("u")
plt.legend()
plt.savefig('plot_betaError.pdf', dpi=300, bbox_inches='tight')
plt.close()