from cProfile import label
from fileinput import filename
import numpy as np
import os
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

def custom_sort(s):
    # Extract the integer part from the string, or use 0 if there are no integers
    integer_part = int(''.join(filter(str.isdigit, s))) if any(char.isdigit() for char in s) else 0
    return (integer_part, s)

dir_list = os.listdir("./")
filtered_dir = [dir_name for dir_name in dir_list if dir_name.startswith('k')]
filtered_dir = sorted(filtered_dir, key=custom_sort) 

colors = ['r', 'b', 'g', 'k']

for i in range(len(filtered_dir)):
    current_directory = os.getcwd()
    current_directory+='/'+filtered_dir[i]
    

    file_list = os.listdir(current_directory)
    

    filtered_files = [filename for filename in file_list if filename.startswith('data_')]
    filtered_files = sorted(filtered_files, key=custom_sort) 

    
    N = np.array([int(filename[5:9]) for filename in filtered_files])
    
    a = 0.0
    b = 1.0
    h = np.array([abs(b-a)/n for n in N])


    L2error = np.genfromtxt(current_directory + '/errorL2.csv', delimiter=',', max_rows=1)

    x = -np.log(h)
    y = np.log(L2error)

    model = LinearRegression()
    x_reshaped = [[val] for val in x]

    h2off = 6

    model.fit(x_reshaped[h2off:], y[h2off:].tolist())
    slope = model.coef_[0]
    slope = abs(slope)
    # print(f"Linear regression slope (angular coefficient): {slope}")

    # alpha = abs(y[-1] - y[-2]) / abs(x[-1] - x[-2])
    # print(alpha)


    plt.plot(x,y,linestyle='-',marker='o', c=colors[i], label='k = ' + filtered_dir[i][1:] + ' | slope = ' + "{:.4f}".format(slope))

    plt.plot(x[h2off:],model.predict(x_reshaped[h2off:]),linestyle='--', c=colors[i])

plt.xlabel("$-\log(h)$")
plt.ylabel("$\log(error)$")
plt.legend()
plt.grid(True)
plt.savefig('ErrorsL2.pdf', dpi=300, bbox_inches='tight')
plt.show()