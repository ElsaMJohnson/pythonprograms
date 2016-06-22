#!/home/elsa/Ureka/variants/common/bin/python
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
import numpy as np
#I imagine this name already has a use, but my program
# is for replotting your final 3d galaxy from the threedgalaxy.txt
#file:
#again designed to cut and paste into ipython --pylab
#read in file - because mixed character and number data can't use genfromtxt
x=[]
y=[]
color=[]
z=[]
fig=figure()
with open('threedgal.txt') as f:
    for line in f:
        line=line.strip()
        col=line.split()
	x.append(float(col[0]))
	y.append(float(col[1]))
        color.append(col[2])	
	z.append(float(col[3]))
	

ax = fig.gca(projection='3d')
nx=max(x)
ny=max(y)
xi = np.linspace(1.0,nx,nx)
yi = np.linspace(1.0,ny,ny)
XX,YY = np.meshgrid(xi,yi)
ZZ = griddata(x,y,ht,xi,yi)
colors = np.array(color).reshape(ny,nx)
surf = ax.plot_surface(XX, YY, ZZ, rstride=1, cstride=1, facecolors=colors,linewidth=0, antialiased=True)
ax.view_init(azim = 270, elev = 90)
ax.set_zlim3d(0,max(ht))
ax.w_zaxis.set_major_locator(LinearLocator(6))
plt.show()

