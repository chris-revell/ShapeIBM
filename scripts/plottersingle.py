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
time = float(argv[5])

grid0 = np.genfromtxt("output/gridpositions0.txt")
grid1 = np.genfromtxt("output/gridpositions1.txt")
data0 = np.genfromtxt("output/fluidvelocities0.txt",skip_header=(Numg*nplot),max_rows=Numg)
data1 = np.genfromtxt("output/fluidvelocities1.txt",skip_header=(Numg*nplot),max_rows=Numg)
data = np.genfromtxt("output/boundarypositions.txt",delimiter=", ",skip_header=(nbounds*nplot),max_rows=nbounds)

fig,ax=plt.subplots()

grad1D = np.linspace(0,0,Numg)
grad2D = np.zeros((Numg,Numg))
for j in range(Numg):
    grad2D[:,j] = grad1D
ax.contourf(grid0,grid1,grad2D,cmap='inferno',levels=Numg)

if colourflag==1:
    M3 = np.sqrt(np.power(data0,2)+np.power(data1,2))
    A = 5*np.divide(data0[::arraystep,::arraystep],M3[::arraystep,::arraystep],out=np.zeros_like(data0[::arraystep,::arraystep]),where=M3[::arraystep,::arraystep]!=0)
    B = 5*np.divide(data1[::arraystep,::arraystep],M3[::arraystep,::arraystep],out=np.zeros_like(data1[::arraystep,::arraystep]),where=M3[::arraystep,::arraystep]!=0)
    ax.quiver(grid0[::arraystep,::arraystep],grid1[::arraystep,::arraystep],A,B,M3[::arraystep,::arraystep],pivot='mid',cmap="Greys")
elif colourflag==0:
    ax.quiver(grid0[::arraystep,::arraystep],grid1[::arraystep,::arraystep],data0[::arraystep,::arraystep],data1[::arraystep,::arraystep],pivot='mid')
ax.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
ax.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)
ax.set_xlim([-5,5])
ax.set_ylim([-5,5])
ax.axis('equal')
ax.plot(np.append(data[:,0],data[0,0]),np.append(data[:,1],data[0,1]))
fig.savefig("output/plot{:04d}.png".format(nplot),bbox_inches='tight',padding_inches=0,dpi=200)

os.system("convert output/plot{:04d}.png -shave 268x143 output/plot{:04d}.png".format(nplot,nplot))
