import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pylab
from numpy import *
from mpl_toolkits.mplot3d import Axes3D
#import cv2.cv as cv
import time
from matplotlib import animation
import subprocess
#plt.use('Agg')


N=100

plot_galaxy= False
galaxy={}
if plot_galaxy: 
	num, posx, posy, posz = np.loadtxt('positions_N100_sm.txt',\
        		usecols=(0,1,2,3), unpack=True) 

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	for i in range(N):
		galaxy[str(i)+'x'] = posx[i::N]
		galaxy[str(i)+'y'] = posy[i::N]
		galaxy[str(i)+'z'] = posz[i::N]
		ax.plot(galaxy[str(i)+'x'], galaxy[str(i)+'y'], galaxy[str(i)+'z'])

	ax.set_ylim([-40,40])
	ax.set_xlim([-40,40])
	ax.set_zlim([-40,40])

	ax.set_title('Open Galactic Cluster in 3D')
	ax.set_xlabel('X-Position [ly]')
	ax.set_ylabel('Y-Position [ly]')
	ax.set_zlabel('Z-Position [ly]')
	plt.show()

#Energy Figures
energy_plot=False
if energy_plot:
	energy, KE, PE = np.loadtxt('Energy_N10.txt',\
        				usecols=(0,1,2), unpack=True) 

	time=np.linspace(1,5,len(energy))

	plt.figure(1)
	plt.plot(time, energy/10000, label='Total Energy')
	plt.plot(time, KE/10000, label='Kinetic Energy')
	plt.plot(time, PE/10000, label='Potential Energy')
	plt.xlabel(r'$\tau_{crunch}$', size=20)
	plt.ylabel('Energy', size=15)
	plt.title('Energy variation with time', size=15)
	plt.legend()
	plt.show()

energy_bounded=False
if energy_bounded:

	PE, KE = np.loadtxt('boundEnergies_Unsmooth200.txt',\
        				usecols=(0,1), unpack=True) 

	time=np.linspace(1,5,len(KE))

	Etot= KE+PE

	plt.figure(1)
	plt.plot(time, KE, label='Kinetic Energy')
	plt.plot(time, Etot, label='Total Energy')
	plt.plot(time, PE, label='Potential Energy')
	plt.xlabel(r'$\tau_{crunch}$', size=20)
	plt.ylabel('Energy', size=14)
	plt.title('Energy Variation With Time N=100 Smoothed', size=16)
	plt.legend()
	plt.show()

#stuff,x,y,z = np.loadtxt('positions_N100_sm', usecols=(0,1,2,3), unpack=True) 

NNbound = True
if NNbound:
	Ntotal = np.array([100.,200.,300.,400.,500.])
	Mass = 1000./Ntotal
	bound = np.array([86.,165.,264.,325.,409.])
	unbound = Ntotal - bound
	fracunbound = unbound/Ntotal
	fracmassejected = (Mass*unbound)/1000.

	print(fracunbound)
	print(mean(fracunbound))

	m, b = np.polyfit(Ntotal, unbound, 1)

	plt.figure(1)
	plt.plot(Ntotal,unbound, 'o', label='data')
	plt.plot(Ntotal, m*Ntotal + b, '-', label='linefit, y=%g x %g' %(m,b))
	plt.xlim(0,550)
	plt.ylim(0,100)
	plt.xlabel(r'Total number of particles N', size=16)
	plt.ylabel('Number of ejected particles', size=16)
	plt.title('Particle ejection plot', size=16)
	plt.legend(loc='best')

	plt.figure(2)
	plt.plot(Ntotal,fracunbound)
	plt.ylim(0,1)
	plt.xlabel(r'Total number of particles N', size=16)
	plt.ylabel('Fraction of ejected particles', size=16)
	plt.title('Fraction of ejected particles', size=16)
	plt.show()

	plt.figure(3)
	plt.plot(Ntotal,fracmassejected)
	plt.ylim(0,1)
	plt.xlabel(r'Total number of particles N', size=16)
	plt.ylabel('Fraction of mass ejected', size=16)
	plt.title('Fraction of mass ejected', size=16)
	plt.show()






