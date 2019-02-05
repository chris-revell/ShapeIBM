# coding: utf-8
import matplotlib.pyplot as plt
import numpy as np
data = np.genfromtxt("output/volume.txt")
fig,ax = plt.subplots()
ax.plot(data[:,0],data[:,1])
ax.set_ylim([0,np.max(data[:,1])])
ax.set_ylabel("Volume /$\mu m^2$")
ax.set_xlabel("Time /s")
fig.savefig("output/volume.png",bbox_inches="tight",padding_inches=0)
