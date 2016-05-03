#!/home/elsa/Ureka/variants/common/bin/python
#Code to generate an image with perlin noise as background. 
#Used in UO scientific computing course Spring 2016
#Perlin code and noise package is from Casey Duncan
#https://github.com/caseman/noise/examples/2dtexture.py 
#Remaining code by Elsa M. Johnson

from noise import pnoise2, snoise2
import numpy as np
from matplotlib import pyplot as plt
import random as random
from PIL import Image
import matplotlib.image as mpimg


#Create noise - 
#size is image size
#octaves creates different sized noise regions
def nbackground(sx=1000,sy=1000,octaves=50):
#	size=1000
#	array = np.zeros((size, size), np.float)
	array = np.zeros((sx, sy), np.float)
#	octaves = 50
	freq = 16.0 * octaves
	for y in xrange(sy):
	    for x in xrange(sx):
        	data=snoise2(x / freq, y / freq, octaves) * 127.0 + 128.0
           	array[x, y] += data


	return array
	plt.imshow(array, cmap=plt.cm.Greys_r)

#To get the pasting right, use cpaste, paste isn't working.

#Now creating a small image note that the image must be smaller than sx and sy
#make sure scspiral.png is in directory
#Will later make the size more variable
#the 
def mkimg(infile='scspiral.png',sz=1000):
	data =[]
#	infile = 'scspiral.png'
	im=Image.open(infile).convert('L')
	plt.imshow(im)
	imlist = list(im.getdata())
	x=int(im.size[1])
	y=im.size[0]
	im.format
	data=np.array(imlist).reshape(x,y)
	cdata=data[50:150,50:150] #cropping image a bit
	#pad with zeros to fit noise image:
	xx,yy=cdata.shape  #size of object
	#randomly pick a beginning number for 
	#location of image and pad with zeros
	xlim = sz-xx-10
	ylim = sz-yy-10
#	These were specific numbers based on the image
#	Which in the end might be better than random placement
#	begx=436
#	begy=596
	begx=random.randint(1,xlim)
	begy=random.randint(1,ylim)
	print 'center of embedded image',begx+50,begy+50 
#	Create limits to pad with zeros
	zx1 = begx-1
	zx2 = sz-zx1-xx
	zy1 = begy-1
	zy2 = sz-zy1-yy
	bimg=np.lib.pad(cdata,((zx1,zx2),(zy1,zy2)),'constant',constant_values=(0,0))
	return bimg


#	This combines both images and scales the added image based on the S/N ratio
#	imarr is the image from mkimg and noisearr is from nbackground
#	sz = is the size of box for doing S/N calculations and s2n is the desired
#	S/N ratio
def combineimg(imarr,noisearr,s2n=5.0,thresh=100,sz=10,):
	b=imarr
	b[b<thresh]=0
	x,y=np.where(b==b.max())
	sig=b[x[0]-5:x[0]+5, y[0]-5:y[0]+5].sum()
	nse=noisearr[x[0]-5:x[0]+5, y[0]-5:y[0]+5].sum()
	#quadratic formula to find fct so that S/N is correct
	fs = (s2n*s2n +np.sqrt(s2n**4.+4*nse*s2n*s2n))/2
	fct=sig/fs
	b=b/fct
	#note to find location of b max: where(b==b.max())
	totimg = b+noisearr
	plt.figure()
	plt.imshow(totimg,cmap=plt.cm.Greys_r)
	return totimg


#Next routine calculates the mean and stand dev of random pixels
def imgstats(arr,sz=100):
#	imm=np.mean(arr[x-sz/2:x+sz/2,y-sz/2:y+sz/2])
#	ims=np.std(arr[x-sz/2:x+sz/2,y-sz/2:y+sz/2])
        ax,ay=arr.shape
	begx=np.random.randint(1,ax-sz)
	begy=np.random.randint(1,ay-sz)
	rm = np.mean(arr[begx:begx+sz,begy:begy+sz])
        rs = np.std(arr[begx:begx+sz,begy:begy+sz])
#	print 'mean,std about image', imm,ims
	print 'random center, mean, std',begx,begy,rm,rs

	#previous values	
	#np.mean(stuff[436:536,596:696])
	#np.std(stuff[436:536,596:696])
	#np.mean(stuff[461:511,621:671])
	#np.std(stuff[461:511,621:671])


def svimg(totarr):
	#print it out:
	x,y=totarr.shape
	vl = np.around(totarr.flatten(),5)#round to 5 digits
	xx = np.repeat(np.arange(x),x)+1
	yy = np.tile(np.arange(y),y)+1
	big =np.column_stack((xx,yy,vl))
	np.savetxt("noisyimage.txt",big,fmt=('%4.1f','%4.1f','%10.5f'))
	##Add this if you want to
	##read it out to make sure it works
	##Otherwise slows down routine.
	#row,col,data=np.loadtxt("noisyimage.txt",unpack=True)
	#rsize = int(max(row))
	#csize = int(max(col))
	#data=np.array(data).reshape(rsize,csize)
#	plt.imshow(data, interpolation='None',cmap=plt.cm.Greys_r)

def main(): 
	noiseimg = nbackground()
	hiddenimg = mkimg()
	timg = combineimg(hiddenimg,noiseimg)
	imgstats(timg)
	svimg(timg)

main()
plt.show()
	
