# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

#### CONSTANTS ####


a = 1.0
b = 1.0



def xpos(a, E, e):
    return a * (np.cos(E) - e)
    
def ypos(b, E, e):
    return b * np.sin(E)
    
Es = np.linspace(0, 2*np.pi, 20)


fig = plt.figure()
ax = fig.add_subplot(111)
plt.scatter(xpos(a, Es, 0.9), ypos(b, Es, 0.9))
plt.scatter(0.,0., color='red')
plt.title("Keplerian Orbit at eccentricity = 0.9")
ax.set_xlabel('X Position')
ax.set_ylabel('Y Position')

plt.savefig("K_e4")

plt.show()