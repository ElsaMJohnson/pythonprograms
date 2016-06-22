#!/usr/bin/python
import numpy as np
import sys, csv

# Just the function, need to make commandline interactive.
# This is to create grids for polyphot
# Box is the size of the area on image
# xcenter ycenter is the center of the box
# res is the pixel size of the grid
# Box needs to contain an integer number of grids and resolution
# overlap is if the grids overlap eachother by 1/2 of the grid
def grids(box=500,xcenter=512,ycenter=512,res=10,overlap='No'):
	xmin= xcenter-box/2+res/2
	ymin = ycenter-box/2+res/2
	xmax = xcenter+1+box/2-res/2
	ymax = ycenter+1+box/2-res/2
	if (overlap=='Yes'):
	    spc = res/2
	elif (overlap =='No'):
	    spc = res
	
	nn = box/spc
	xx = np.tile(np.arange(xmin,xmax,spc),nn)
	yy = np.repeat(np.arange(ymin,ymax,spc),nn)
	big =np.column_stack((xx,yy))
	np.savetxt("grids.txt",big,fmt=('%4.0f','%4.0f'))

