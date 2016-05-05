import numpy as np

nx = 10 #number of grids in x and y
ny = 10 
xbeg = 0 #x and y start pixels
ybeg = 0 
xend=101  #x and y end pixels
yend = 101
spc = 10  #space between grid centers
xx = np.tile(np.arange(xbeg,xend,spc),ny)+1
yy = np.repeat(np.arange(ybeg,yend,spc),nx)+1
big =np.column_stack((xx,yy))
np.savetxt("grids.txt",big,fmt=('%4.0f','%4.0f'))

