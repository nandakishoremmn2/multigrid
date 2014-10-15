'''
	Plots residue vs work units data
	Usage:
		python plot_data.py [input file]

	'input fie' default = 'data.dat'
'''

import sys

from pylab import *

data = loadtxt('data.dat' if len(sys.argv) < 2 else sys.argv[1])

plot(cumsum(data[:,2]), data[:,1], '-o')
yscale('log')
axis('image')

show()