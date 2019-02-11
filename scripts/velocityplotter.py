import numpy as np
import matplotlib.pyplot as plt
import os

Numg=512
colourflag=1

os.system("rm output/velocity*.png;rm output/velocityanimated.gif;")

grid0 = np.genfromtxt("output/gridpositions0.txt")
grid1 = np.genfromtxt("output/gridpositions1.txt")
data0 = np.genfromtxt("output/fluidvelocities0.txt")
data1 = np.genfromtxt("output/fluidvelocities1.txt")
data = np.genfromtxt("output/boundarypositions.txt",delimiter=", ")
nbounds = np.genfromtxt("output/nbounds.txt",dtype=int)

fig,ax=plt.subplots()
nsteps = int(np.shape(data0)[0]/(Numg+1))
drawn = 0
for i in range(nsteps):
    fig,ax=plt.subplots()
    if colourflag==1:
        M3 = np.sqrt(np.power(data0,2)+np.power(data1,2))
        A = 5*np.divide(data0[i*(Numg+1):(i+1)*(Numg+1):8,::8],M3[i*(Numg+1):(i+1)*(Numg+1):8,::8],out=np.zeros_like(data0[i*(Numg+1):(i+1)*(Numg+1):8,::8]), where=M3[i*(Numg+1):(i+1)*(Numg+1):8,::8]!=0)
        B = 5*np.divide(data1[i*(Numg+1):(i+1)*(Numg+1):8,::8],M3[i*(Numg+1):(i+1)*(Numg+1):8,::8],out=np.zeros_like(data1[i*(Numg+1):(i+1)*(Numg+1):8,::8]), where=M3[i*(Numg+1):(i+1)*(Numg+1):8,::8]!=0)
        ax.quiver(grid0[::8,::8],grid1[::8,::8],A,B,M3[i*(Numg+1):(i+1)*(Numg+1):8,::8],pivot='mid',cmap="Greys")
    else:
        ax.quiver(grid0[::8,::8],grid1[::8,::8],data0[i*(Numg+1):(i+1)*(Numg+1):8,::8],data1[i*(Numg+1):(i+1)*(Numg+1):8,::8],pivot='mid')
    ax.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
    ax.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)
    ax.set_xlim([-5,5])
    ax.set_ylim([-5,5])
    ax.axis('equal')
    ax.plot(np.append(data[drawn:drawn+nbounds[i],0],data[drawn,0]),np.append(data[drawn:drawn+nbounds[i],1],data[drawn,1]))
    fig.savefig("output/velocitytest{:04d}.png".format(i),bbox_inches='tight',padding_inches=0,dpi=200)
    plt.close()
    drawn = drawn+nbounds[i]
os.system("for i in $(ls output|grep velocity); do convert output/$i -shave 268x143 output/$i; done;")
os.system("convert -delay 10 -loop 0 output/velocitytest*.png output/velocityanimated.gif;")
