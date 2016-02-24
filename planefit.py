import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


def read_fits(datafile):

    import pyfits
    import os
 
    if not os.path.exists(datafile):
        print "***error! data file",data_file," does not exist!"
        return 0
    hdulist = pyfits.open(datafile)
    return np.asarray(hdulist[1].data)

   
    
      
datafile = 'C:\Users\Administrator\Desktop\Astro\Data\karachentsev13.csv'
data = np.genfromtxt(datafile, names=True, delimiter=',', dtype=None)
glon = data['ra']
glat = data['dec']
#bmag = data['Vmag']#ABS_BMAG
D = data['distance'] # heliocentric distance in kpc




glen = len(glon)
glong = [0.0] * glen
glatt = [0.0] * glen
DD = [0.0] * glen
mag = [0.0] * glen


for i in range(glen - 1):
    glong[0]=0
    glatt[0]=0
    DD[0]=0
    glong[i] = float(glon[i+1])
    glatt[i] = float(glat[i+1])
    DD[i+1] = float(D[i+1])/1000# * 1000 convert to kpc if necessary (Heasarc)
    mag[i] = 1#float(0.02 * np.power((np.sqrt(np.power(bmag[i],2))),np.e))
        
gl = [0.0] * glen
gt = [0.0] * glen
SGB = [0.0] * glen
K = [0.0] * glen
Q = [0.0] * glen
J =  [0.0] * glen
SGL =  [0.0] * glen
SGx =  [0.0] * glen
SGy =  [0.0] * glen
SGz =  [0.0] * glen
sgc =  [0.0] * glen
    
def SGconvert(glon,glat,D): #converts from galactic to supergalactic 3D
    SGNra = np.radians(47.37)
    SGNdec = np.radians(6.32)
    SGCra = np.radians(137.37)
    SGCdec = 0.0
    for i in range (len(glon)):
        gl[i] = ((glong[i]) / 360) * 2 * np.pi
        gt[i] = ((glatt[i]) / 360) * 2 * np.pi
        SGB[i] = np.arcsin(np.sin(gt[i]) * np.sin(SGNdec) 
        + np.cos(gt[i]) * np.cos(SGNdec) * np.cos(gl[i] - SGNra))
        K[i] = np.arcsin((np.cos(gt[i]) * np.sin(gl[i] - SGNra)) 
        / np.cos(SGB[i]))
        J[i] = (np.sin(gt[i]) * np.cos(SGNdec) 
        - np.cos(gt[i])*np.sin(SGNdec)*np.cos(gl[i] - SGNra)) / np.cos(SGB[i])
        Q[i] = np.arccos(np.sin(SGCdec) / np.cos(SGNdec))
        
        if J[i] < 0:
            SGL[i] = Q[i] + K[i] - np.radians(180)
        else:
                SGL[i] = Q[i] - K[i]
                
        #if SGL[i] < 0:
         #   SGL[i] += np.radians(360)
            
        SGy[i] = DD[i] * np.cos(SGB[i]) * np.sin(SGL[i])
        SGx[i] = DD[i] * np.cos(SGB[i]) * np.cos(SGL[i])
        SGz[i] = DD[i] * np.sin(SGB[i])
                            
    sgc = zip(SGx,SGy,SGz)
    return np.asarray(sgc)

##Plot SG co-ordinates
SGconvert(glong,glatt,DD)
from mayavi import mlab
mlab.outline(mlab.axes(mlab.points3d(SGx, SGy, SGz)))
mlab.show()

#plt.scatter(SGx,SGy)
fig = plt.scatter(SGz,SGy,s = mag)
plt.ylabel('SGy / Mpc')
plt.xlabel('SGz / Mpc')
plt.xlim((-11,11))
plt.ylim((-11,11))
plt.show()


def fitPlaneSVD(XYZ):
    [rows,cols] = XYZ.shape
    # Set up constraint equations of the form  AB = 0,
    # where B is a column vector of the plane coefficients
    # in the form b(1)*X + b(2)*Y +b(3)*Z + b(4) = 0.
    p = (np.ones((rows,1)))
    AB = np.hstack([XYZ,p])
    [u, d, v] = np.linalg.svd(AB,0)        
    B = v[3,:];                    # Solution is last column of v.
    nn = np.linalg.norm(B[0:3])
    B = B / nn
    return B
    

   
#fitPlaneSVD(SGconvert(glong,glatt,DD))
 

b = fitPlaneSVD(SGconvert(glong,glatt,DD))
print (SGconvert(glong,glatt,DD))
print b


def Drms(sgc, b):
    dist = [0] * glen
    Dsq = 0
    sgx = sgc[:,0]
    sgy = sgc[:,1]
    sgz = sgc[:,2]
    for i in range(glen):
        dist[i] = (b[0]*sgx[i] + b[1]*sgy[i] + b[2]*sgz[i] + b[3]) / np.sqrt(np.square(b[0]) + np.square(b[1]) + np.square(b[2]))
        Dsq += np.square(dist[i])
    return np.sqrt(Dsq / glen)

print " "
print "The Drms for this sample = %.2f Mpc" % (Drms(SGconvert(glong,glatt,DD), b))
print " "
 
print " "
print "fitting ``fundamental'' plane to the sample"
print "    log(Re)=a*mu + b*log(sigma) + const"
print 'plane fit results: a = %.2f, b = %.2f' % (-b[1]/b[0], -b[2]/b[0])
print " "


#print "NOW TESTING"
#testco = [(0,0,0)] * 10
#testco[9] = (10,-10,10)
#for i in range (1,10):
#    testco[i-1] = (i * np.power(-1,i),i * np.power(-1,i),i * np.power(-1,i))
#testc = np.asarray(testco)
#tb = [0]*4
#tb[2] = 1 #tb is a plane on z=0
#print testc
#print tb
#    
#    
#def testDrms(testc, b):
#    dist = [0] * 10
#    Dsq = 0
#    testx = testc[:,0]
#    testy = testc[:,1]
#    testz = testc[:,2]
#    print b[2]
#    print testz
#    for i in range(10):
#        dist[i] = (b[0]*testx[i] + b[1]*testy[i] + b[2]*testz[i] + b[3]) / np.sqrt(np.square(b[0]) + np.square(b[1]) + np.square(b[2]))
#        print dist[i]
#        Dsq += np.square(dist[i])
#    return np.sqrt(Dsq / 10)
#    
#print " "
#print "The Drms for TEST sample = %.2f" % (testDrms(testc,tb))
#print " "


#bootstrapping
import numpy as np
from numpy import random as rnd

x = SGconvert(glong,glatt,DD)
#
# now bootstrap
#
nboot = 10000

from scipy.optimize import curve_fit

def flin(x,m,c):
    return m*x + c

mbchain = []

#mm = 0; sm = 0;
#for ib in range(nboot):
#    xb = np.zeros_like(x)
#    for i in range(len(x)):
#        ii = rnd.randint(len(x)-1)
#        xb[i]=x[ii]
#    b0 = fitPlaneSVD(xb)[0]
#    b1 = fitPlaneSVD(xb)[1]
#    b2 = fitPlaneSVD(xb)[2]
#    #mb = zip(b0,b1,b2)
#    mm += mb; sm += mb**2
#    mbchain.append(mb)
#    
#print "bootstrap results:"
#mm = mm/nboot; sm = np.sqrt(sm/nboot-mm**2)
#print "m e_m c e_c=",mm, sm

#sb2 = np.sum((x-flin(x,mm,sm))**2)/len(x) - mm**2*np.sum(ex**2)/len(x)
#print "scatter", np.sqrt(sb2)

# results of a simple least squares fit

#p, cov = curve_fit(fitPlaneSVD,x,sigma=ex)
#mlsq = p[0]; clsq = p[1]
#print "results of a simple least squares fit"
#print "m, em, c, ec=", mlsq, np.sqrt(cov[0,0]), clsq, np.sqrt(cov[1,1])
