import numpy as np
import matplotlib.pyplot as plt
import os
from sys import argv

nplot = int(argv[1])
nbounds = int(argv[2])
Numg=int(argv[3])+1
colourflag=int(argv[4])

grid0 = np.genfromtxt("output/gridpositions0.txt")
grid1 = np.genfromtxt("output/gridpositions1.txt")
data0 = np.genfromtxt("output/fluidvelocities0.txt")
data1 = np.genfromtxt("output/fluidvelocities1.txt")
data = np.genfromtxt("output/boundarypositions.txt",delimiter=", ")

fig,ax=plt.subplots()

if colourflag==1:
    M3 = np.sqrt(np.power(data0,2)+np.power(data1,2))
    A = 5*np.divide(data0[-Numg::8,::8],M3[-Numg::8,::8],out=np.zeros_like(data0[-Numg::8,::8]),where=M3[-Numg::8,::8]!=0)
    B = 5*np.divide(data1[-Numg::8,::8],M3[-Numg::8,::8],out=np.zeros_like(data1[-Numg::8,::8]),where=M3[-Numg::8,::8]!=0)
    ax.quiver(grid0[-Numg::8,::8],grid1[-Numg::8,::8],A,B,M3[-Numg::8,::8],pivot='mid',cmap="Greys")
else:
    ax.quiver(grid0[::8,::8],grid1[::8,::8],data0[-Numg::8,::8],data1[-Numg::8,::8],pivot='mid')
ax.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
ax.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)
ax.set_xlim([-5,5])
ax.set_ylim([-5,5])
ax.axis('equal')
ax.plot(np.append(data[-nbounds:,0],data[-nbounds,0]),np.append(data[-nbounds:,1],data[-nbounds,1]))
fig.savefig("output/velocitytest{:04d}.png".format(nplot),bbox_inches='tight',padding_inches=0,dpi=200)

os.system("convert output/velocitytest{:04d}.png -shave 268x143 output/velocitytest{:04d}.png".format(nplot,nplot))
