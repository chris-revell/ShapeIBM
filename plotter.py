import numpy as np
import matplotlib.pyplot as plt
import os

data = np.genfromtxt("boundarypositions.txt",delimiter=", ")

xmax=10
xmin=10
npoints = 64
nsteps = int(np.shape(data)[0]/npoints)
for step in range(nsteps):
    fig, ax = plt.subplots(figsize=(4,4))
    ax.scatter(data[step*npoints:(step+1)*npoints,0],data[step*npoints:(step+1)*npoints,1])
    ax.set_xlim([-5,5])
    ax.set_ylim([-5,5])
    ax.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
    ax.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)
    fig.savefig("test{:02d}".format(step))
    plt.close()
os.system("convert -delay 10 -loop 0 *.png animated.gif; rm *.png")
