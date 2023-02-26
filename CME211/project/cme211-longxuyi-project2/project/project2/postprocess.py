import sys
import math
import numpy as np
import warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import interp1d

# print usage message if missing arguments
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print('Usage:')
        print('  python3 postprocess.py <input_file> <soln_file>')
        sys.exit(0)

input_file = sys.argv[1]
soln_file = sys.argv[2]

# read input file
f = open(input_file, 'r')
size = f.readline()
l, w, h = float(size.split()[0]), float(size.split()[1]), float(size.split()[2])
temp = f.readline()
Tc, Th = float(temp.split()[0]), float(temp.split()[1])
f.close()

# create a 2d temperature array
nrow = int (w/h + 1)
ncol = int (l/h + 1)
temp_data = np.zeros([nrow, ncol])

#add cold boundary into the 2d array
for j in range(ncol):
    temp_data[0][j] = -Tc*(math.exp(-10*(j*h - l/2)*(j*h - l/2)) - 2)

#read solution file into the temperature matrix
f = open(soln_file, 'r')
for i in range(1,nrow-1):
    for j in range(ncol):
        if j == 0:
            periodic_temp = f.readline().split()[0]
            temp_data[i][j] = periodic_temp
        elif j == ncol-1:
            temp_data[i][j] = periodic_temp
        else:
             temp_data[i][j] = f.readline().split()[0]
f.close()

#add hot boundary into the 2d array
for j in range(ncol):
    temp_data[nrow-1][j] = Th

temp_data = np.flipud(temp_data)
mean =np.mean(temp_data[:, 1:])
print("Input file processed: {}".format(input_file))
print("Mean Temperature: {}".format('%.5f' % mean))


# plot heat map and color bar
plt.imshow(temp_data, extent = (0.0, l, 0.0, w))
plt.colorbar()

# plot mean temp isoline
x = np.linspace(0,l,ncol)
y = np.linspace(0,w,nrow)
plt.contour(x, w - y, temp_data - mean, levels = [0], linewidth = 2, colors = "black")

plt.title("Pipe Wall Temperature Distribution")
plt.xlabel('x')
plt.ylabel('y')
plt.xlim(0, l)
plt.ylim(0, w)
plt.show()

plt.savefig('test')