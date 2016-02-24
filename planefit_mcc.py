import numpy as np
from matplotlib import pylab as plt

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
 
def read_McConnachie12():
    datafile = 'C:\Users\Administrator\Desktop\Astro\Data\mcconnachie12.fits'
    data = read_fits(datafile)
    # ra and dec 
    RA = data['RAJ2000']; DEC = data['DEJ2000']
    # heliocentric distance in kpc
    D = data['D']
    d = D.astype(int)
    d[0] = 0; d = d.astype(float)/1000. # kpc-> Mpc
    ra = np.zeros_like(d); dec = np.zeros_like(d)
    for i in range(len(RA)):
        ras = np.fromstring(RA[i], sep=" ")
        decs = np.fromstring(DEC[i], sep=" ")
        ra[i] = (ras[0]*3600.*15.+ ras[1]*60. + ras[2])*np.pi/(180.*3600) # 360/24 = 15
        dec[i] = np.sign(decs[0])*(np.abs(decs[0])*3600.+ decs[1]*60. + decs[2])*np.pi/(180.*3600)
        #print "%3d %2.1f %2.1f %2.1f %2.4f %2.4f %4.4f"%(i, ras[0], ras[1], ras[2], ra[i], dec[i], d[i])
    return ra, dec, d
    
def read_Karachentsev13():
    datafile = 'C:\Users\Administrator\Desktop\Astro\Data\karachentsev13.fits'
    data = read_fits(datafile)
    # ra and dec 
    RA = data['RA']; DEC = data['DEC']
    # heliocentric distance in kpc
    D = data['DISTANCE'] 
    #VPEC = data['VPEC']
    d = D.astype(float)
    ras = RA.astype(float); decs = DEC.astype(float)
    ra = ras * np.pi/180.; dec = decs * np.pi/180.
    return ra, dec, d
    
def plot_pretty():
    plt.rc('text', usetex=True)
    plt.rc('font',size=20)
    plt.rc('xtick.major',pad=5); plt.rc('xtick.minor',pad=5)
    plt.rc('ytick.major',pad=5); plt.rc('ytick.minor',pad=5)

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
    
def plot_SG_projections(SGXk, SGYk, SGZk, SGXm, SGYm, SGZm):
    plot_pretty()
 
    plt.figure(figsize=(17.0,8.0))
    plt.subplot(1,2,1)
    plt.xlabel(r'$\rm SGZ\ (Mpc)$'); plt.ylabel(r'$\rm SGY\ (Mpc)$')
    plt.xlim(-11.,11.); plt.ylim(-11.,11.)
    plt.scatter(SGZk, SGYk)
    plt.scatter(SGZm, SGYm, c='r')

    plt.grid()
    plt.show()

    plt.subplot(1,2,2)
    plt.xlabel(r'$\rm SGX\ (Mpc)$'); plt.ylabel(r'$\rm SGY\ (Mpc)$')
    plt.xlim(-11.,11.); plt.ylim(-11.,11.)
    plt.scatter(SGXk, SGYk, c='b', label=r'$\rm Karachentsev\ et\ al.\ 2012$')
    plt.scatter(SGXm, SGYm, c='r', label=r'$\rm McConnachie\ et\ al.\ 2012$')

    plt.grid()
    plt.legend(frameon=False, loc='lower right', fontsize=14)
    plt.show()
    
ram, decm, dm = read_McConnachie12()

#isel = (dm < 0.250)
#ram = ram[isel]; decm =decm[isel]; dm = dm[isel]
#ram = ram[1:]; decm = decm[1:]; dm=dm[1:]

#rak, deck, dk = read_Karachentsev13()
#isel = (dk<11.5)
#rak = rak[isel]; deck =deck[isel]; dk = dk[isel]
#rak = rak[1:]; deck = deck[1:]; dk=dk[1:]
#xk, yk, zk = eq_to_3d(rak, deck)
xm, ym, zm = eq_to_3d(ram, decm)


SGXm, SGYm, SGZm = dm*eq2000_to_supergal(xm, ym, zm)

"""
SGXk, SGYk, SGZk = dk*eq2000_to_supergal(xk,yk,zk)


XYZk = np.array(zip(SGXk, SGYk, SGZk))
bk = fitPlaneSVD(XYZk)
 
print " "
print "fitting plane to the K13 sample"
print "    a*SGX + b*SGY + c*SGZ + d=0"
print 'plane fit results: a = %.2f, b = %.2f, c=%.2f, d=%.2f' % (bk[0], bk[1], bk[2], bk[3])
print " "
"""

XYZm = np.array(zip(SGXm, SGYm, SGZm))
XYZmb = np.array(zip(SGXm, SGYm, SGZm))
bm = fitPlaneSVD(XYZm)
 
print " "
print "fitting plane to the McConnachie12 sample"
print "    a*SGX + b*SGY + c*SGZ + d=0"
print 'plane fit results: a = %.2f, b = %.2f, c=%.2f, d=%.2f' % (bm[0], bm[1], bm[2], bm[3])
print " "

dpm = (bm[0]*SGXm + bm[1]*SGYm + bm[2]*SGZm + bm[3])/np.sqrt(bm[0]**2 + bm[1]**2 + bm[2]**2)

Drmsm = np.std(dpm)

print "Drmsm = ", Drmsm


"""
# calculate angle between the planes defined for two samples
nk = bk[:3]; nm = bm[:3]

angle = np.arccos(np.dot(nk,nm)/np.sqrt(nk.dot(nk))/np.sqrt(nm.dot(nm))) * 180./np.pi
print "angle between the two planes is ", angle, " degrees"


# check out distribution of galaxies in projections visually
plot_SG_projections(SGXk, SGYk, SGZk, SGXm, SGYm, SGZm)

"""


#___BOOTSTRAPPING___
#McConnachie
nboot = 10000
dist = [0] * len(SGXm)
bb = [0] * nboot
Drmsmb = [0] * nboot
angle = [0] * nboot

def mcbootstrap(SGXYZ,Drms,B):
    Drmsum = 0
    for i in range(nboot):
        Dsq = 0
        for k in range(len(SGXm)):
            ii = np.random.randint(len(SGXm)-1)
            XYZmb[k] = SGXYZ[ii]
        bb[i] = fitPlaneSVD(XYZmb)
        for j in range(len(SGXm)):
            dist[j] = (bb[i][0]*XYZmb[j][0] + bb[i][1]*XYZmb[j][1] + bb[i][2]*XYZmb[j][2] + bb[i][3])/np.sqrt(bb[i][0]**2 + bb[i][1]**2 + bb[i][2]**2)
            Dsq += np.square(dist[j])
        Drmsmb[i] = np.sqrt(Dsq / len(SGXm))
        nk = B[:3]
        nm = bb[i][:3]
        angle[i] = np.arccos(np.dot(nk,nm)/np.sqrt(nk.dot(nk))/np.sqrt(nm.dot(nm))) * 180./np.pi
        if angle[i] > 90:
            angle[i] = -(180 - angle[i])
        Drmsum += np.square(Drmsmb[i] - Drms)
    Drmstdv = np.sqrt(Drmsum / nboot)
    anglestdv = np.std(angle)
    return Drmstdv, anglestdv
        
Drmstd, anglestd = mcbootstrap(XYZm,0.41288554488,fitPlaneSVD(XYZm))
print "Drms stdev =", Drmstd
print "angle stdev =", anglestd

