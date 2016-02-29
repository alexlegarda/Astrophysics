import numpy as np
import matplotlib.pyplot as plt

###READ DATA###

def read_fits(data_file=None):
    """Loader for a fits table of galaxies
    
    Returns
    -------
    data : recarray, shape = (327260,)
        record array containing pipeline parameters
 
 
    """
    # pyfits is an optional dependency: don't import globally
    import pyfits
    import os
 
    if not os.path.exists(data_file):
        print "***error! data file",data_file," does not exist!"
        return 0
    hdulist = pyfits.open(data_file)
    return np.asarray(hdulist[1].data)

def read_McConnachie12():
    datafile = r'data/mcconnachie12.fits'
    data = read_fits(datafile)
    mcgalname = data['Name']
    Vmag = data['VMag1']
    vmag = []
    for i in range(len(Vmag)):
        vmag.extend(np.fromstring(Vmag[i], sep=" "))
    #from operator import itemgetter
    #print min(enumerate(vmag), key=itemgetter(1))[0]
    #print mcgalname[0]
    # ra and dec 
    RA = data['RAJ2000']; DEC = data['DEJ2000']
    # heliocentric distance in kpc
    D = data['D_MW_']
    d = D.astype(int)
    d[0] = 0; d = d.astype(float)/1000. # kpc-> Mpc
    
    # V-band luminosity
    #for i in range(len(vmag)):
        #vmag[i] = vmag[i] + 5 - 5*np.log10(d[i])
    lv = [0]*len(vmag)
    for i in range(len(vmag)):
        # solar luminosity
        lv[i] = pow(10, (5.48 - vmag[i])/2.5)
        
    vr = data['V_MW_']; vr[0] = 0.;
    ra = np.zeros_like(d); dec = np.zeros_like(d)
    for i in range(len(RA)):
        ras = np.fromstring(RA[i], sep=" ")
        decs = np.fromstring(DEC[i], sep=" ")
        if vr[i] == '    ':
            vr[i] = 0.
        else:
            vr[i] = vr[i].astype(float)
        ra[i] = (ras[0]*3600.*15.+ ras[1]*60. + ras[2])*np.pi/(180.*3600) # 360/24 = 15
        dec[i] = np.sign(decs[0])*(np.abs(decs[0])*3600.+ 
        decs[1]*60. + decs[2])*np.pi/(180.*3600)
        #print "%3d %s %2.1f %2.1f"%(i, mcgalname[i], d[i], vr[i])
        #print i, mcgalname[i], d[i], vr[i]
    vr = vr.astype(float)
    print 'mcconachie sample size =', len(ra)
    return ra, dec, d, lv
    
def read_Karachentsev13():
    datafile = r'data\karachentsev13_BMag.txt'
    data = np.loadtxt(datafile, skiprows=1, delimiter=',')
    # ra and dec 
    RA = data[:,0]; DEC = data[:,1]
    # heliocentric distance in kpc
    D = data[:,2]
    #VPEC = data['VPEC']
    #Bolometric Magnitude
    bm = data[:,3]
    L = [0]*len(bm)
    m = [0]*len(bm)
    #Mass calculation
    for i in range(len(bm)):
        # solar luminosity
        L[i] = pow(10, (5.48 - bm[i])/2.5) # http://www.ucolick.org/~cnaw/sun.html
        #if L[i] > 5e10:
            #print D[i], bm[i], L[i]
    d = D.astype(float)
    ras = RA.astype(float); decs = DEC.astype(float)
    ra = ras * np.pi/180.; dec = decs * np.pi/180.
    
    #D Selection#
    dmax = 30 #MODIFY SELECTION HERE, MAX is 26.2
    Ra, Dec, Dd, M, Ll = zip(*((ra, dec, d, m, L) 
    for ra, dec, d, m, L in zip(ra, dec, d, m, L) if d<dmax))
    print 'Karachentsev D Selection:', max(Dd), '<', dmax, 'Mpc'  #Sanity check
    return Ra, Dec, Dd, Ll

def eq2000_to_supergal(x, y, z):
    """
    convert from equatorial to supergalactic coordinates
    using rotation matrix multiplication. 
    for rotation matrices fror different transformation see:
        http://www.atnf.csiro.au/computing/software/gipsy/sub/skyco.c
    input: vectors of x,y,z coordinates obtained from equatorial ra, dec as
        x = np.cos(ra) * np.cos(dec)
        y = np.sin(ra) * np.cos(dec)
        z = np.sin(dec)
    return: s = array of triples of transformed coordinates
        multiply them by distance to get SGX, SGY, SGZ
    
    """
    Rot = [
      [0.3751891698,   0.3408758302,   0.8619957978],
      [-0.8982988298,  -0.0957026824,   0.4288358766],
      [0.2286750954,  -0.9352243929,   0.2703017493]]     
      
    r = np.array([x,y,z])
    s = np.zeros_like(r)
    s = np.dot(Rot,r)
    return s
    
def eq1950_to_supergal(x, y, z):
    Rot = [
      [0.3831583954,   0.3366379840,   0.8601537722],
      [-0.8972185056,  -0.0856688522,   0.4331971849],
      [0.2195190133,  -0.9377290203,   0.2692130889]]
    r = np.array([x,y,z])
    s = np.zeros_like(r)
    s = np.dot(Rot,r)
    return s
    
def eq_to_3d(ra, dec):
    x = np.cos(ra) * np.cos(dec)
    y = np.sin(ra) * np.cos(dec)
    z = np.sin(dec)
    return x, y, z



## READ AND MERGE DATA ##
ramm, decmm, dmm, Lmm = read_McConnachie12()
ram, decm, dm, Lm = read_Karachentsev13()

xm, ym, zm = eq_to_3d(ram, decm)
xmm, ymm, zmm = eq_to_3d(ramm, decmm)
SGXm, SGYm, SGZm = dm*eq2000_to_supergal(xm, ym, zm)
SGXmm, SGYmm, SGZmm = dmm*eq2000_to_supergal(xmm, ymm, zmm)
SGS = zip(SGXm,SGYm,SGZm,Lm)
print 'Karachentsev dataset length =', len(SGS)
SGSm = zip(SGXmm,SGYmm,SGZmm,Lmm)
seen = set(item[:3] for item in SGS)
SGS.extend(item for item in SGSm[:3] if item[:3] not in seen)
print 'total dataset length =', len(SGS)
print 'data added from mcconachie:', SGS[-3:]


## FINAL DATA ARRAYS ##
SGXm, SGYm, SGZm, m = np.asarray(zip(*SGS))
m[-3] = pow(10, (5.48 - 19.1)/2.5) #use Karachentsev magnitude for MW

#SGXm -= np.average(SGXm); SGYm -= np.average(SGYm); SGZm -= np.average(SGZm)
ic = [0.,0.,0.] #search center
Rm = np.sqrt((SGXm - ic[0])**2 + (SGYm - ic[1])**2 + (SGZm - ic[2])**2)
Rmax=5.
#Select by distance from search center
SGXm, SGYm, SGZm, m, Rm = np.asarray(zip(*((SGXm, SGYm, SGZm, m, Rm) 
    for SGXm, SGYm, SGZm, m, Rm in zip(SGXm, SGYm, SGZm, m, Rm) if Rm<Rmax)))
CofMX, CofMY, CofMZ = 0, 0, 0
for i in range(len(SGXm)):
    CofMX += SGXm[i]*m[i]
    CofMY += SGYm[i]*m[i] 
    CofMZ += SGZm[i]*m[i]
CofMX /= np.sum(m); CofMY /= np.sum(m); CofMZ /= np.sum(m)
print 'SG Radial Selecion at center [', ic[0], ',', ic[1], ',', ic[2], ']' 
print 'and radius:', max(Rm), '<', Rmax, 'Mpc'
print 'SG CofM of sample = [',CofMX,',', CofMY,',', CofMZ,']'
print 'sum of masses:', np.sum(m)
print 'max m:', max(m), 'min m:', min(m), 'mean m:', np.mean(m)


# test the calculation with synthetic data points
"""
sigx = 3.0; sigy = 3.0; sigz = 0.03
npt = 1000
SGXm = np.random.normal(0.0,sigx,npt)
SGYm = np.random.normal(0.0,sigy,npt)
SGZm = np.random.normal(0.0,sigz,npt)
"""



XYZm = np.array(zip(SGXm, SGYm, SGZm))
#**TEST**#
###REDEFINE XYZ and m FOR TESTING PURPOSES###

#MANUAL
#xs = (0, 1, -1, 0, 0, 0, 0)
#ys = (0, 0, 0, 1, -1, 0, 0)
#zs = (0, 0, 0, 0, 0, 1, -1)
#XYZm = np.array(zip(xs, ys, zs))
#m = np.array([1.0] * len(XYZm))
#
##RANDOM
#tl = 100
#XYZm = 10*np.random.rand(tl, 3)
#m = np.random.rand(len(XYZm), 1)[:,0]
##m[3] = 10
#
#CofMX = 0.
#CofMY = 0.
#CofMZ = 0.
#for i in range(len(XYZm)):
#    print XYZm[:,0][i]
#    print m[i]
#    CofMX += XYZm[:,0][i]*m[i]
#    CofMY += XYZm[:,1][i]*m[i] 
#    CofMZ += XYZm[:,2][i]*m[i]
#CofMX /= np.sum(m); CofMY /= np.sum(m); CofMZ /= np.sum(m)
#print 'SG CofM of TEST sample = [',CofMX,',', CofMY,',', CofMZ,']'


###ELLIPSOID FIT###



def ellfit(XYZ, m):#EIGENVECTORS
    from numpy import linalg as LA
    nxx = 0.0
    nyy = 0.0
    nzz = 0.0
    nxy = 0.0
    nyz = 0.0
    nzx = 0.0
    x = XYZ[:,0]
    y = XYZ[:,1]
    z = XYZ[:,2]
    M = np.sum(m)
    for i in range(len(XYZ)):
       nxx = np.sum(m * x * x) / M #nxx += m[i] * x[i] * x[i] / M
       nyy += m[i] * y[i] * y[i] / M
       nzz += m[i] * z[i] * z[i] / M
       nxy += m[i] * x[i] * y[i] / M
       nyz += m[i] * y[i] * z[i] / M
       nzx += m[i] * x[i] * z[i] / M
    S = np.matrix([[nxx,nxy,nzx],[nxy,nyy, nyz],[nzx,nyz,nzz]])
    eigenValues, eigenVectors = LA.eig(S)
    return eigenVectors, eigenValues
    
def ellfit2(XYZ):#EIGENVALUES
    from numpy import linalg as LA
    #m = [1.0] * len(XYZ) #TO BE MODIFIED
    #for i in range(len(XYZ)): 
      #  m[i] = i
    r = [0] * len(XYZ)
    x = XYZ[:,0]
    y = XYZ[:,1]
    z = XYZ[:,2] 
    M = np.sum(m)
    S = (3,3); S = np.zeros(S)
    for i in range(len(XYZ)):
        S[0,0] += m[i] * x[i] * x[i] / M
        S[1,1] += m[i] * y[i] * y[i] / M
        S[2,2] += m[i] * z[i] * z[i] / M
        S[0,1] += m[i] * x[i] * y[i] / M
        S[1,2] += m[i] * y[i] * z[i] / M
        S[0,2] += m[i] * x[i] * z[i] / M
    S[1,0] = S[0,1]; S[2,1] = S[1,2]; S[2,0]=S[0,2]
    eigenValues = LA.eigvals(S)
    return eigenValues
    
### BOOTSTRAPPING ###

boot = np.array(zip(XYZm[:,0],XYZm[:,1],XYZm[:,2],m))
it = 10000
bootvals = [0] * it
bootvecs = [0] * it
for i in range(it):
    strap = np.array(boot[np.random.randint(len(boot), size = len(boot)),:])
    bootvecs[i], bootvals[i] = ellfit(np.array(zip(strap[:,0], strap[:,1], strap[:,2])), np.array(strap[:,3]))
bootvals = np.array(bootvals)

xyboot = bootvals[:,0] / bootvals[:,1]
zxboot = bootvals[:,2] / bootvals[:,0]
zyboot = bootvals[:,2] / bootvals[:,1]
xystd = np.std(xyboot)
zxstd = np.std(zxboot)
zystd = np.std(zyboot)
xymean = np.mean(xyboot)
zxmean = np.mean(zxboot)
zymean = np.mean(zyboot)

print ' '
print 'Bootstrapping results:'
print 'iterations =', len(bootvals)
print 'given as (mean, stdev)'
print 'X/Y---------------Z/X---------------Z/Y'
print('(%.4f,%.4f),(%.4f,%.4f),(%.4f,%.4f)'%
(xymean,xystd,zxmean,zxstd,zymean,zystd))
print ' '




#m = np.array([1.0] * len(XYZm))
eigenVectors, eigenValues = ellfit(XYZm,m)
sqrteV = np.sqrt(eigenValues)
print "sqrt eigenvalues =", sqrteV
#eigenValues = ellfit2(XYZm)
print "eigenvalues =", eigenValues

#ellispsoid and center in matrix form
A = np.array(eigenVectors)
center = [CofMX,CofMY,CofMZ]
print "eigenvectors=\n", eigenVectors

radii = 3*sqrteV

print "ellipsoid radii:", radii
print "Z/X =", (sqrteV[2]/sqrteV[0])
print "Z/Y =", (sqrteV[2]/sqrteV[1])
print "X/Y =", (sqrteV[0]/sqrteV[1])
rotation = A.T



# assign ellipsoid parameters
u = np.linspace(0.0, 2.0 * np.pi, 100)
v = np.linspace(0.0, np.pi, 100)
x = radii[0] * np.outer(np.cos(u), np.sin(v))
y = radii[1] * np.outer(np.sin(u), np.sin(v))
z = radii[2] * np.outer(np.ones_like(u), np.cos(v))
for i in range(len(x)):
    for j in range(len(x)):
        [x[i,j],y[i,j],z[i,j]] = np.dot([x[i,j],y[i,j],z[i,j]],rotation)+center

# plot
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()

ax = Axes3D(fig)
#rlim = max(max(SGXm), max(SGYm), max(SGZm))
#ax.set_xlim(-rlim,rlim); ax.set_ylim(-rlim,rlim); ax.set_zlim(-rlim,rlim)
ax.set_zlabel('SGZ')
ax.set_xlabel('SGX')
ax.set_ylabel('SGY')

ax.plot_wireframe(x, y, z,  rstride=4, cstride=4, color='r', alpha=0.2)

plt.hold
print "max(m)=", max(m)
s = m/max(m)*500
ax.scatter(SGXm,SGYm,SGZm, s=s)
#ax.scatter(XYZm[:,0], XYZm[:,1], XYZm[:,2], s=s)

###Bounding box + limits to set equal aspect ratio for axes
max_range = np.array([SGXm.max()-SGXm.min(), SGYm.max()-SGYm.min(), SGZm.max()-SGZm.min()]).max() / 2.0
mean_x = SGXm.mean()
mean_y = SGYm.mean()
mean_z = SGZm.mean()
ax.set_xlim(mean_x - max_range, mean_x + max_range)
ax.set_ylim(mean_y - max_range, mean_y + max_range)
ax.set_zlim(mean_z - max_range, mean_z + max_range)

plt.show()
#plt.close(fig)
#del fig
