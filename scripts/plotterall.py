import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
from sys import argv
from math import sqrt

Numg=int(argv[1])+1
arraystep=int(Numg/64)
xmin = -1 #float(argv[3])
xmax = 1 #float(argv[4])

grid0 = np.genfromtxt(argv[2]+"/gridpositions0.txt")
grid1 = np.genfromtxt(argv[2]+"/gridpositions1.txt")
nbounddata = np.genfromtxt(argv[2]+"/nbounds.txt")
data = np.genfromtxt(argv[2]+"/boundarypositions.txt",delimiter=", ")

nboundspassed = 0

for nplot,floatnbounds in enumerate(nbounddata):
    nbounds = int(floatnbounds)
    fig,ax=plt.subplots()
    ax.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
    ax.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)

    ordered = np.zeros((nbounds,2))
    i=0
    currentelement = 0
    previouselement = 0
    while i<nbounds:
        neighbour = int(data[nboundspassed+currentelement,-2])
        if neighbour == previouselement:
            neighbour = int(data[nboundspassed+currentelement,-1])
        ordered[i,:] = data[nboundspassed+neighbour,:2]
        i=i+1
        previouselement=currentelement
        currentelement = neighbour


    ax.plot(np.append(ordered[:,0],ordered[0,0]),np.append(ordered[:,1],ordered[0,1]),c="black")
    #ax.scatter(data[nboundspassed:nboundspassed+nbounds,0],data[nboundspassed:nboundspassed+nbounds,1],c=data[nboundspassed:nboundspassed+nbounds,2],s=5,cmap='gray_r',vmin=0,vmax=np.max(data[:,2]),edgecolors='none')
    ax.set_xlim([xmin,xmax])
    ax.set_ylim([xmin,xmax])
    ax.set_aspect('equal')
    nboundspassed = nboundspassed+nbounds
    fig.savefig(argv[2]+"/plot{:04d}.png".format(nplot),bbox_inches='tight',padding_inches=0,dpi=200)
    plt.close()

#os.system("convert "+argv[2]+"/plot*.png -shave 268x143 "+argv[2]+"/plot*.png")
