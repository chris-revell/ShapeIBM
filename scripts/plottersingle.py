import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
from sys import argv
from math import sqrt

nplot = int(argv[1])
nbounds = int(argv[2])
Numg=int(argv[3])+1
colourflag=int(argv[4])
arraystep=int(Numg/64)
xmin = float(argv[5])
xmax = float(argv[6])
time = float(argv[7])

grid0 = np.genfromtxt(argv[8]+"/gridpositions0.txt")
grid1 = np.genfromtxt(argv[8]+"/gridpositions1.txt")
data0 = np.genfromtxt(argv[8]+"/fluidvelocities0.txt",skip_header=(Numg*nplot),max_rows=Numg)
data1 = np.genfromtxt(argv[8]+"/fluidvelocities1.txt",skip_header=(Numg*nplot),max_rows=Numg)
data = np.genfromtxt(argv[8]+"/boundarypositions.txt",delimiter=", ",skip_header=(nbounds*nplot),max_rows=nbounds)

fig,ax=plt.subplots()

if colourflag==1:
    M3 = np.sqrt(np.power(data0,2)+np.power(data1,2))
    A = 5*np.divide(data0[::arraystep,::arraystep],M3[::arraystep,::arraystep],out=np.zeros_like(data0[::arraystep,::arraystep]),where=M3[::arraystep,::arraystep]!=0)
    B = 5*np.divide(data1[::arraystep,::arraystep],M3[::arraystep,::arraystep],out=np.zeros_like(data1[::arraystep,::arraystep]),where=M3[::arraystep,::arraystep]!=0)
    ax.quiver(grid0[::arraystep,::arraystep],grid1[::arraystep,::arraystep],A,B,M3[::arraystep,::arraystep],pivot='mid',cmap="Greys")
elif colourflag==0:
    ax.quiver(grid0[::arraystep,::arraystep],grid1[::arraystep,::arraystep],data0[::arraystep,::arraystep],data1[::arraystep,::arraystep],pivot='mid')
ax.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
ax.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)
#ax.plot(np.append(data[:,0],data[0,0]),np.append(data[:,1],data[0,1]))
ax.scatter(data[:,0],data[:,1],c=data[:,2],s=5,cmap='gray_r',vmin=0,vmax=np.max(data[:,2]),edgecolors='none')
ax.set_xlim([xmin,xmax])
ax.set_ylim([xmin,xmax])
ax.set_aspect('equal')
fig.savefig(argv[8]+"/plot{:04d}.png".format(nplot),bbox_inches='tight',padding_inches=0,dpi=200)

#os.system("convert output/plot{:04d}.png -shave 268x143 output/plot{:04d}.png".format(nplot,nplot))
