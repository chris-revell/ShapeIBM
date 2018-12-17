import numpy as np
import matplotlib.pyplot as plt
import os

os.system("rm output/*.png;rm output/volanimated.gif;")

grid0 = np.genfromtxt("output/gridpositions0.txt")
grid1 = np.genfromtxt("output/gridpositions1.txt")
data0 = np.genfromtxt("output/fluidvelocities0.txt")
data1 = np.genfromtxt("output/fluidvelocities1.txt")
data = np.genfromtxt("output/boundarypositions.txt",delimiter=", ")
nbounds = np.genfromtxt("output/nbounds.txt",delimiter=", ",dtype=int)

fig,ax=plt.subplots()

nsteps = int(np.shape(data0)[0]/65)
drawn = 0
for i in range(nsteps):
    if i%10<1:
        ax.quiver(grid0[::2,::2],grid1[::2,::2],data0[i*65:(i+1)*65:2,::2],data1[i*65:(i+1)*65:2,::2])
        ax.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
        ax.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)
        ax.plot(np.append(data[drawn:drawn+nbounds[i],0],data[drawn,0]),np.append(data[drawn:drawn+nbounds[i],1],data[drawn,1]))
        fig.savefig("output/voltest{:04d}.png".format(i),bbox_inches='tight',padding_inches=0,dpi=200)
        ax.cla()
    drawn = drawn+nbounds[i]
os.system("convert -delay 10 -loop 0 output/voltest*.png output/volanimated.gif;")
