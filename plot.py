'''
	Generates surface plot of input file
	Usage:
		python plot.py [input file]

	'input fie' default = 'out.dat'
'''

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from numpy import loadtxt, linspace, meshgrid
import sys

fig = plt.figure()
# ax = fig.add_subplot(111)
ax = fig.add_subplot(111, projection='3d')

z = loadtxt('out.dat' if len(sys.argv) < 2 else sys.argv[1])

n = z.shape[0]

x, y = meshgrid(linspace(0,1,n), linspace(0,1,n))

ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='b')
# ax.contourf(x, y, z,  100, rstride=4, cstride=4, color='b')
ax.axis('image')

plt.show()

