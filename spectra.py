#Reading SDSS spectra (put data into a table and plot)
#Will have to come up with a fitter to flatten data.
#data is in m39/spectra
import pyfits
import math
import numpy
from matplotlib import pylab as plt
file="spec-2240-53823-0566.fits"
fl = pyfits.open(file)
data = fl[1].data
cols=fl[1].columns
nrows=data.shape
flux = data['flux']
ll=data['loglam']
wl = 10**ll
ax=plt.axes()
ax.plot(wl,flux,'-k')
plt.show()

