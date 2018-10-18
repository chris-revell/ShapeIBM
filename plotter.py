import numpy as np
import matplotlib.pyplot as plt


data = np.genfromtxt("boundarypositions.txt",delimiter=", ")

npoints = 64
nsteps = int(np.shape(data)[0]/npoints)
for step in range(nsteps):
    fig, ax = plt.subplots(figsize=(4,4))
    ax.scatter(data[step*npoints:(step+1)*npoints,0],data[step*npoints:(step+1)*npoints,1])
    plt.show()
