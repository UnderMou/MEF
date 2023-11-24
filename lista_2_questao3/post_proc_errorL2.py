from fileinput import filename
import numpy as np
import os
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

def custom_sort(s):
    # Extract the integer part from the string, or use 0 if there are no integers
    integer_part = int(''.join(filter(str.isdigit, s))) if any(char.isdigit() for char in s) else 0
    return (integer_part, s)


file_list = os.listdir("./")

filtered_files = [filename for filename in file_list if filename.startswith('data_')]
filtered_files = sorted(filtered_files, key=custom_sort) 
N = np.array([int(filename[5:9]) for filename in filtered_files])

a = 0.0
b = 1.0
h = np.array([abs(b-a)/n for n in N])


L2error = np.genfromtxt('errorL2.csv', delimiter=',', max_rows=1)
print(N,h)
print(L2error)

x = -np.log(h)
y = np.log(L2error)

# model = LinearRegression()
# x_reshaped = [[val] for val in x]
# model.fit(x_reshaped[0:], y[0:].tolist())
# slope = model.coef_[0]
# slope = abs(slope)
# print(f"Linear regression slope (angular coefficient): {slope}")

alpha = abs(y[-1] - y[-2]) / abs(x[-1] - x[-2])
print(alpha)


plt.plot(x,y,linestyle='-',marker='o')
plt.xlabel("$-\log(h)$")
plt.ylabel("$\log(error)$")
plt.grid(True)
plt.show()