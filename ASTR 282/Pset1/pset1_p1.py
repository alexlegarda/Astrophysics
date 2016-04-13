from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt

def meanM(E,e):
    return E - e * np.sin(E)
    
N = 100

x = np.linspace(0,20,N)

e1 = [0.0] * N
e2 = [0.2] * N
e3 = [0.6] * N
e4 = [0.9] * N

fig = plt.figure()
ax = fig.add_subplot(111)
plt.scatter(meanM(x, e4), x)
plt.title("Eccentric Anomaly vs. Mean Anomaly at eccentricity = 0.9")
ax.set_xlabel('Mean Anomaly')
ax.set_ylabel('Eccentric Anomaly')

plt.savefig("e4")

plt.show()