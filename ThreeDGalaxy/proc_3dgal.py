#!/home/elsa/Ureka/variants/common/bin/python
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
import numpy as np

##By Elsa Johnson
##Heat filter algorithm by Sam Gerber (originally written for MATLAB)
##
##following MATLAB functions did not exist and had to be written in python:
##MATLAB function  ---> Python function
##fspecial --> fspecialgauss2D  
## --from stackoverflow http://stackoverflow.com/questions/17190649/how-to-obtain-a-gaussian-filter-in-python
##fspecial --> fspecialdisk

##Note, routine has many parts designed to just cut and paste into ipython environment.
## Also start ipython with the --pylab flag. Just makes it easier.

##This version take the output put from readpolyphot and then calculates br and puts
## it into the br steps, plots various stages, save to file and then shows the 3d version!
# args:
# 	file - input file
# 	res - size of each grid box/data point in units of pixels
#	cut - color cut off. for example, anything smaller than the cutoff will be considered blue and plotted as blue and anything
#	     bigger than cutoff will be plotted red. I try to set this to Sun colors.
#	Blim and rlim are the faintest surface brightnesses of each grid and anything fainter in either band will be 'green'
#       overlap is in reference to the grids. If yes, that means the grids from polyphot overlap by 1/2 a grid
#	outfile is outfile
##Example usage
#Z,Xn,Yn,color=subgrids('m81medreadpoly.txt',res=25,cut=.6, blim=23,rlim=22,overlap='No', outfile = 'unfilteredM81.txt')
#heatfilter(Z,val=2,outfile='m81filtered.txt')


def subgrids(file, res=25, cut=1, blim=23, rlim=22,overlap = 'No',outfile='out.txt'):
#  galfile='m101output24nol.txt'
#res=24.
   if (overlap == 'Yes'):
	area = 4*res*res
   elif (overlap == 'No'):
	area = res*res
   sbfactor = 2.5*np.log10(area)
   xx,yy,bm,rm,ber,rer = np.genfromtxt(file,dtype = float, skiprows=1,unpack=True)
   Xn=np.divide(np.subtract(xx,min(xx)),res) +1.
   Yn=np.divide(np.subtract(yy,min(yy)),res) +1
   mX=max(Xn)
   mY=max(Yn)
   xlen = len(Xn)
   ylen = len(Yn)
   br=[]
   color = []
   Z = []
   b = bm+sbfactor
   r = rm+sbfactor
   br=np.subtract(b,r)
#Thinking about the cut off. If red is really red b/c blue is really faint
# Need to change this as per resolution:
# for 30 res (15, this was 23.1 and 22 resp)
   for i in range(len(br)):
         if br[i]<=cut:
                if b[i]>=blim:
                        color.append('g')
                        Z.append(0.)
                elif b[i]<blim:
                        color.append('b')
                        Z.append(abs(round(b[i],1)-blim))
         elif br[i]>cut:
                if r[i]>=rlim:
                        color.append('g')
                        Z.append(0.)
                elif r[i]<rlim:
                        color.append('r')
                        Z.append(abs(round(r[i],1)-rlim))

#if you want to save to a text file at this point:
   np.savetxt(outfile,np.column_stack((Xn,Yn,Z)),fmt=('%5.2f','%5.2f','%10.5f'))
   Z = np.array(Z).reshape(mX,mY)
   plt.figure()
   #Below is for color
   #imshow(Z)
   #This is for grayscale
   plt.imshow(Z,cmap = cm.Greys_r)
   return Z,Xn,Yn,color

##Heat filter to get rid of stars
##First all of the special subroutines
import math
import scipy.ndimage as nd
from skimage.morphology import disk
def fspecialgauss2D(shape=(3,3),sigma=0.5):
    """
    2D gaussian mask - should give the same result as MATLAB's
    fspecial('gaussian',[shape],[sigma])
    """
    m,n = [(ss-1.)/2. for ss in shape]
    y,x = np.ogrid[-m:m+1,-n:n+1]
    h = np.exp( -(x*x + y*y) / (2.*sigma*sigma) )
 ## The statement below determines the machine
 ## epsilon - if gaussian is smaller than that
 ## set to 0.
    h[ h < np.finfo(h.dtype).eps*h.max() ] = 0
    sumh = h.sum()
 ## This notation states that it takes the input
 ## h and then it divides by the sum and returns h 
    if sumh != 0:
        h /= sumh
    return h

def fspecialdisk(r=3.0):
    """
    This is matlab code translated into python code for fspecialdisk:
    """
    cr = math.ceil(r-0.5)
    x,y = np.mgrid[-cr:cr+1,-cr:cr+1]
    mxxy = np.maximum(abs(x),abs(y))
    mnxy = np.minimum(abs(x),abs(y))
    m1 = np.zeros((2*r+1,2*r+1))
    m1 = (r**2.< (mxxy+.5)**2. +(mnxy-.5)**2.)*(mnxy-0.5)+ \
     np.nan_to_num((r**2.>=(mxxy+.5)**2. + \
     (mnxy-.5)**2.)*sqrt(r**2.-(mxxy+0.5)**2.))
    m2 = (r**2.>(mxxy-.5)**2. +(mnxy+.5)**2.)*(mnxy+0.5)+ \
     (r**2.<=(mxxy-.5)**2. +(mnxy+.5)**2.)*np.sqrt(r**2.-(mxxy-0.5)**2.)
    sgrid = ((r**2.)*(0.5*(arcsin(m2/r)-arcsin(m1/r))+ \
     0.25*(sin(2*arcsin(m2/r))-sin(2*arcsin(m1/r))))- \
     (mxxy-0.5)*(m2-m1)+(m1-mnxy+0.5))\
     *(logical_or(logical_and((r**2.<(mxxy+.5)**2.+(mnxy+0.5)**2.), \
     (r**2.>(mxxy-0.5)**2.+ (mnxy-0.5)**2.)),\
     (((mnxy==0)&(mxxy-0.5<r)&(mxxy+.5>=r))))).astype(float)
    sgrid =sgrid+((mxxy+.5)**2.+(mnxy+0.5)**2.<r**2.)
    sgrid[cr,cr] = min(pi*r**2,pi/2)
    if cr>0 and r>cr-0.5 and r**2.< (cr-0.5)**2.+0.25:
        m1 =sqrt(r**2-(cr-.5)**2.)
        m1n=m1/r
        sg0=2*(r**2.*(0.5*arcsin(m1n)+0.25*sin(2*arcsin(m1n)))-m1*(cr-0.5))
        sgrid[2*cr,cr]=sg0
        sgrid[cr,2*cr]=sg0
        sgrid[cr,0]=sg0
        sgrid[0,cr]=sg0
        sgrid[2*cr-1,cr]=sgrid[2*cr-1,cr]-sg0
        sgrid[cr,2*cr-1]=sgrid[cr,2*cr-1]-sg0
        sgrid[cr,1]=sgrid[cr,1]-sg0
        sgrid[1,cr]=sgrid[1,cr]-sg0
        
    sgrid[cr,cr]=min(sgrid[cr,cr],1)
    h=sgrid/sgrid.sum()
    return h


def heatfilter(Z,val=.6,outfile='filterout.txt' ):
## Here, image is Z data
   im=Z
   G=fspecialgauss2D((3,3),0.7)
##Apply gaussian filter
   imB = nd.correlate(im,G,mode='reflect')
##Operate laplacian (del2 in matlab)
   L=nd.filters.laplace(imB)
##Set all negative 2nd derivatives to 0
#This was commented out for some reason
#L(L<0)=0
##Use new gaussian filter
   G=fspecialgauss2D((3,3),1)
##Apply this filter
   L=nd.correlate(abs(L),G,mode='reflect')
##Call it a mask
   mask=L
#mask[mask<15]=0
   mask[mask<val]=0
   mask[mask>0] = 1;
#Create a structuring element of a disk;
   sel = disk(2)
#Apply structuring element to mask with imdilate
   mask = nd.morphology.binary_dilation(mask,sel).astype(mask.dtype)
   X=im
   X[mask>0]=0
#X=(X/X.max())*255.##for movie
   xd,yd = X.shape
#movie=zeros((xd,yd,5001)) 
   G = fspecialdisk(1) 
   iter = 0
   delta =1.
   while iter<5001:
        Xg = nd.correlate(X,G,mode='reflect')
        Xtmp = X
        Xtmp[mask>0]=Xg[mask>0]
#        Xtmp = (Xtmp/X.max())*255.
        delta - sum((X[mask>0]-Xtmp[mask>0])**2.)
#        movie[:,:,iter]=X
        iter=iter+1
        X=Xtmp      

#Uncomment this and all other movie references if you want to see your results as a movie
#Takes a lot of memory, so I don't use it.
#movie[:,:,iter:]= []
#movie=movie/255
#Note fspecialgauss2D can be changed: make the arguments larger for bigger grids.
   G=fspecialgauss2D((3,3),0.5)
   X = nd.correlate(X,G,mode='reflect')
   figure()
#for default color image uncomment below
  #imshow(X)
# for gray scale
   imshow(X,cmap = cm.Greys_r)
#Okay now put this back into a file. This file will need to be translated into 3d printerfile
   crap = X.ravel()
# Save to a basic txt file of image with data values:
   savetxt(outfile,column_stack((Xn,Yn,crap)),fmt=('%5.2f','%5.2f','%10.5f'))
   of = open('threedgal.txt','w')
   for i in range(len(crap)):
        of.write(str(Xn[i])+"  "+str(Yn[i])+"  "+str(color[i])+"  "+str(round(crap[i],1))+"\n")

   fig = figure()
   ax = fig.gca(projection='3d')
   nx = max(Xn)-min(Xn) +1.
   ny = max(Yn) - min(Yn) +1.
   xi = np.linspace(min(Xn),max(Xn),nx)
   yi = np.linspace(min(Yn),max(Yn),ny)
   XX,YY = np.meshgrid(xi,yi)
   ZZ = griddata(Xn,Yn,crap,xi,yi)
   colors = np.array(color).reshape(max(Yn),max(Xn))
   surf = ax.plot_surface(XX, YY, ZZ, rstride=1, cstride=1, facecolors=colors,linewidth=0, antialiased=True)
#This is to view from the top #change as necessary
   ax.view_init(azim = 270, elev = 90)
   ax.set_zlim3d(0,max(crap))
   ax.w_zaxis.set_major_locator(LinearLocator(6))
   plt.show()


