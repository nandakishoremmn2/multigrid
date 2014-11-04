'''
	Plots residue vs work units data
	Usage:
		python plot_data.py [input file]

	'input fie' default = 'data.dat'
'''

import sys

from pylab import *

data = loadtxt('data.dat' if len(sys.argv) < 2 else sys.argv[1])


wd = data[:,2]
wu = cumsum(data[:,2])
res = data[:,1]
mask = (data[:,3] == 0)
plot(wu, res, '-o', wu[mask], res[mask], '-o')
yscale('log')
axis('image')
title('Plot of residue')
xlabel('Work Units')
ylabel('norm-2 of residue')
legend(['residue', 'residue of fine grid alone'])

show()

mask = mask * (wd != 0)

wdm = wd[mask]
wum = wu[mask]
resm = res[mask]
plot(wum[1:], pow(resm[1:]/resm[:-1], 1/wdm[1:]), '-o')
title('Smoothing Factor of residue (avgeraged)')
xlabel('Work Units')
ylabel('Smoothing Factor')

show()