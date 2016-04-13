import numpy as np
import matplotlib.pyplot as plt

fs = [300, 99, 49, 15, 2, 3, 1]
ps = [2,3,4,5,6,7,8]
h = [0]*7
for i in range(len(ps)):
    h[i] = [ps[i]]*fs[i]
print h


fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist(h, ps)
plt.title("Planets Per Exoplanet System")
plt.xlabel("Frequency")
plt.ylabel("Number of Planets")
plt.show()
plt.savefig('planets_histogram')