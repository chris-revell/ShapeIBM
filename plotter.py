import numpy as np
import matplotlib.pyplot as plt
import os

data = np.genfromtxt("output/boundarypositions.txt",delimiter=", ")
Nc = 3
xmax=10
npoints = 64*Nc
nsteps = int(np.shape(data)[0]/(npoints))
for step in range(nsteps):
    print("{:02d}/{:02d}".format((step+1),nsteps))
    fig, ax = plt.subplots(figsize=(4,4))
    ax.scatter(data[step*npoints:(step+1)*npoints,0],data[step*npoints:(step+1)*npoints,1])
    ax.set_xlim([-xmax/4,xmax/4])
    ax.set_ylim([-xmax/4,xmax/4])
    ax.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
    ax.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)
    fig.savefig("output/test{:02d}".format(step))
    plt.close()
os.system("convert -delay 10 -loop 0 output/*.png output/animated.gif;")
# rm output/*.png")
