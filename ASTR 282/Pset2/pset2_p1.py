import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize


datafile1 = r'assignment2_prob1_data1.txt'
data1 = np.loadtxt(datafile1, skiprows=2)
t1 = data1[:,0]
f1 = data1[:,1]#/np.max(data1[:,1])
u1 = data1[:,2]#/np.max(data1[:,1])

datafile2 = r'assignment2_prob1_data2.txt'
data2 = np.loadtxt(datafile2, skiprows=2)
t2 = data2[:,0]
f2 = data2[:,1]#/np.max(data2[:,1])
u2 = data2[:,2]#/np.max(data2[:,1])

datafile3 = r'assignment2_prob1_data3.txt'
data3 = np.loadtxt(datafile3, skiprows=2)
t3 = data3[:,0]
f3 = data3[:,1]#/np.max(data3[:,1])
u3 = data3[:,2]#/np.max(data3[:,1])



#Constants

e = 0.0
K = 141.24
Ms = 0.9


def sin(p, x):
    return p[0] * np.sin(x * p[1] + p[2]) + p[3]
def residual(p, x, y, err):
    return (sin(p, x) - y) / err

p1 = [1000.,1.,0.,0.]    
pf1, cov1, info1, mesg1, success1 = optimize.leastsq(residual, p1, args=(t1, f1, u1), full_output=1)
chisq1 = sum(info1["fvec"]*info1["fvec"])
dof1 = len(t1)-len(pf1)
pferr1 = [np.sqrt(cov1[i,i]) for i in range(len(pf1))]

p2 = [1000.,1.,0.,0.]    
pf2, cov2, info2, mesg2, success2 = optimize.leastsq(residual, p2, args=(t2, f2, u2), full_output=1)
chisq2 = sum(info2["fvec"]*info2["fvec"])
dof2 = len(t2)-len(pf2)
pferr2 = [np.sqrt(cov2[i,i]) for i in range(len(pf2))]

p3 = [1000.,10.,0.,0.]    
pf3, cov3, info3, mesg3, success3 = optimize.leastsq(residual, p3, args=(t3, f3, u3), full_output=1)
chisq3 = sum(info3["fvec"]*info3["fvec"])
dof3 = len(t3)-len(pf3)
pferr3 = [np.sqrt(cov3[i,i]) for i in range(len(pf3))]


fig1, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1)

ax1.errorbar(t1, f1, yerr=u1, fmt='b.', label = 'Data')
T = np.linspace(t1.min(), t1.max(), 5000)
#ax1.plot(T, sin(pf1, T), 'r-', label = 'curve')
ax1.set_title('Transit Data Set 1')
ax1.set_xlabel('Time (days)')
ax1.set_ylabel('Flux (normalized)')
textfit1 = 'Period = %.3f$\pm %.3f$ days \n' \
           % (np.pi/pf1[1], (pferr1[1]/pf1[1])*np.pi/pf1[1])
ax1.text(0.02, .3, textfit1, transform=ax1.transAxes, fontsize=12,
         verticalalignment='top')

ax2.errorbar(t2, f2, yerr=u2, fmt='b.', label = 'Data')
T = np.linspace(t2.min(), t2.max(), 5000)
#ax2.plot(T, sin(pf2, T), 'r-', label = 'curve')
ax2.set_title('Transit Data Set 2')
ax2.set_xlabel('Time (days)')
ax2.set_ylabel('Flux (normalized)')
textfit2 = 'Period = %.3f$\pm %.3f$ days \n' \
           % (np.pi/pf2[1], (pferr2[1]/pf2[1])*np.pi/pf2[1])
ax2.text(0.02, .3, textfit2, transform=ax2.transAxes, fontsize=12,
         verticalalignment='top')

ax3.errorbar(t3, f3, yerr=u3, fmt='b.', label = 'Data')
T = np.linspace(t3.min(), t3.max(), 5000)
#ax3.plot(T, sin(pf3, T), 'r-', label = 'curve')
ax3.set_title('Transit Data Set 3')
ax3.set_xlabel('Time (days)')
ax3.set_ylabel('Flux (normalized)')
textfit3 = 'Period = %.3f$\pm %.3f$ days \n' \
           % (np.pi/pf3[1], (pferr3[1]/pf3[1])*np.pi/pf3[1])
ax3.text(0.65, .3, textfit3, transform=ax3.transAxes, fontsize=12,
         verticalalignment='top')

plt.tight_layout()
plt.show()
plt.savefig('transits.png')


Rratio1 = np.sqrt(np.abs(2*pf1[0])/np.max(data1[:,1]))
Rratio1err = Rratio1 * pferr1[0]/pf1[0]
Rratio2 = np.sqrt(np.abs(pf2[0])/np.max(data2[:,1]))
Rratio2err = Rratio2 * pferr2[0]/pf2[0]
Rratio3 = np.sqrt(np.abs(pf3[0])/np.max(data3[:,1]))
Rratio3err = Rratio3 * pferr3[0]/pf3[0]

Rratio = (Rratio1 + Rratio2 + Rratio3) / 3
Rratioerr = Rratio * np.sqrt((Rratio1err/Rratio1)*(Rratio1err/Rratio1)+(Rratio2err/Rratio2)*(Rratio2err/Rratio2)+(Rratio3err/Rratio3)*(Rratio3err/Rratio3))

print Rratio
print Rratioerr

#0.103165876155
#0.00429143513873
