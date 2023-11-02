import numpy as np
import os
import matplotlib.pyplot as plt

def f(x):
    # return np.multiply(np.sin(np.pi*x),np.sin(np.pi*x))
    return np.sin(np.pi*x)

def custom_sort(s):
    # Extract the integer part from the string, or use 0 if there are no integers
    integer_part = int(''.join(filter(str.isdigit, s))) if any(char.isdigit() for char in s) else 0
    return (integer_part, s)

x = np.linspace(-2,2,1024)
u = f(x)

file_list = os.listdir("./")

filtered_files = [filename for filename in file_list if filename.startswith('data_')]
filtered_files = sorted(filtered_files, key=custom_sort) 

for file in filtered_files:

    x_projL2 = np.genfromtxt(file, delimiter=',', max_rows=1)
    u_projL2 = np.genfromtxt(file, delimiter=',', skip_header=1, max_rows=1)

    plt.figure()
    plt.grid(True)

    plt.plot(x,u,c = "b",label="$u(x) = sin(\pi x)$")
    plt.scatter(x_projL2,u_projL2,s=5,c = "r",label="Projeção $L_2$", zorder=2)
    
    
    plt.xlabel("x")
    plt.ylabel("u")
    plt.legend()
    plt.savefig('plot'+file[:-4]+'.pdf', dpi=300, bbox_inches='tight')
    plt.close()

N = np.array([4,8,16,32,64,128,256])
N = np.divide(4,N)
L2error = np.array([1.08884, 0.219947, 0.0525804, 0.0129298, 0.00321801, 0.000803316, 0.000200908])

plt.plot(-np.log(N),np.log(L2error),linestyle='-',marker='o')
plt.xlabel("$$")
plt.show()