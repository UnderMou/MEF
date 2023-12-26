import numpy as np
import os
import matplotlib.pyplot as plt
import math

def f(x, epsilon, kappa):
    fu = np.zeros_like(x)
    for i in range(len(x)):
        fu[i] = -1.90476*math.exp(-5.0*x[i])+3.9604*math.exp(-x[i])+(5.99498*math.exp(100*x[i]))/pow(10,44)-2.05563
    return fu

def custom_sort(s):
    # Extract the integer part from the string, or use 0 if there are no integers
    integer_part = int(''.join(filter(str.isdigit, s))) if any(char.isdigit() for char in s) else 0
    return (integer_part, s)

x = np.linspace(0.0,1.0,1024)
epsilon = 1e-2
kappa = 1.0
u = f(x, epsilon, kappa)

file_list = os.listdir("./")

filtered_files = [filename for filename in file_list if filename.startswith('data_')]
filtered_files = sorted(filtered_files, key=custom_sort) 

for file in filtered_files:

    x_projL2 = np.genfromtxt(file, delimiter=',', max_rows=1)
    u_projL2 = np.genfromtxt(file, delimiter=',', skip_header=1, max_rows=1)

    plt.figure()
    plt.grid(True)

    plt.plot(x,u,c = "b",label="$u(x) exata$")
    plt.plot(x_projL2, u_projL2, c="r", marker='o', linestyle='--', label="MEF", zorder=2)
    
    plt.xlabel("x")
    plt.ylabel("u")
    plt.legend()
    plt.savefig('plot'+file[:-4]+'.pdf', dpi=300, bbox_inches='tight')
    plt.close()

# N = np.array([4,8,16,32,64,128,256])
# N = np.divide(4,N)
# L2error = np.array([1.08884, 0.219947, 0.0525804, 0.0129298, 0.00321801, 0.000803316, 0.000200908])

# plt.plot(-np.log(N),np.log(L2error),linestyle='-',marker='o')
# plt.xlabel("$$")
# plt.show()