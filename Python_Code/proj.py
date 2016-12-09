import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pylab
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


%matplotlib inline
#IMPORTING DATA FROM QT
data = loadtxt("positions_Tau2.txt")
dataJ = loadtxt("positions_Tau2.txt", skiprows=2)
dataM = loadtxt("positions_Tau2.txt", skiprows=3)
dataV = loadtxt("positions_Tau2.txt", skiprows=4)
dataS = loadtxt("positions_Tau2.txt", skiprows=5)
dataMer = loadtxt("positions_Tau2.txt", skiprows=6)
dataU = loadtxt("positions_Tau2.txt", skiprows=7)
dataN = loadtxt("positions_Tau2.txt", skiprows=8)
dataP = loadtxt("positions_Tau2.txt", skiprows=9)
dataposEarth = data[1::10]
dataposE= dataposEarth[1::10]
dataposJupiter = dataJ[::10]
dataposMars = dataM[::10]
dataposVenus = dataV[::10]
dataposSaturn = dataS[::10]
dataposMercury = dataMer[::10]
dataposUranus = dataU[::10]
dataposNeptune = dataN[::10]
dataposPluto = dataP[::10]
xE,yE,zE = dataposE[:,1], dataposE[:,2]  , dataposE[:,3]
xJ,yJ,zJ = dataposJupiter[:,1], dataposJupiter[:,2]  , dataposJupiter[:,3]
xM,yM,zM = dataposMars[:,1], dataposMars[:,2]  , dataposMars[:,3]
xV,yV,zV = dataposVenus[:,1], dataposVenus[:,2]  , dataposVenus[:,3]
xS,yS,zS = dataposSaturn[:,1], dataposSaturn[:,2]  , dataposSaturn[:,3]
xMer,yMer,zMer = dataposMercury[:,1], dataposMercury[:,2]  , dataposMercury[:,3]
xU,yU,zU = dataposUranus[:,1], dataposUranus[:,2]  , dataposUranus[:,3]
xN,yN,zN = dataposNeptune[:,1], dataposNeptune[:,2]  , dataposNeptune[:,3]
xP,yP,zP = dataposPluto[:,1], dataposPluto[:,2]  , dataposPluto[:,3]
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(xE, yE, zE)
ax.plot(xJ, yJ, zJ)
ax.plot(xM, yM, zM)
ax.plot(xV,yV,zV)
ax.plot(xS,yS,zS)
ax.plot(xMer,yMer,zMer)
ax.plot(xU,yU,zU)
ax.plot(xN,yN,zN)
ax.plot(xP,yP,zP)
ax.set_title('Galaxy Cluster in 3D')
ax.set_xlabel('X-Position [AU]')
ax.set_ylabel('Y-Position [AU]')
ax.set_zlabel('Z-Position [AU]')
#ax.scatter([0],[0],[0],color="r",s=150)
#ax.scatter([0.8],[0.6],[0],color="g",s=50)
#ax.scatter([5.1],[0],[0],color="y",s=70)
#ax.set_aspect('equal')
#ax.set_ylim([-10,10])
#ax.set_xlim([-10,10])
#ax.set_zlim([-10,10])
#ax.annotate('Earth',xy=(0.01,0.003))
#ax.annotate('SUN',xy=(-0.002,0.007))
#ax.annotate('Jupiter',xy=(0.03,0.0))
plt.show()