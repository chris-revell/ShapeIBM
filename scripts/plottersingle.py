import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
from sys import argv
from math import sqrt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

nplot = int(argv[1])
nbounds = int(argv[2])
Numg=int(argv[3])+1
colourflag=int(argv[4])
arraystep=int(Numg/64)
xmin = float(argv[5])
xmax = float(argv[6])

grid0 = np.genfromtxt(argv[7]+"/gridpositions0.txt")
grid1 = np.genfromtxt(argv[7]+"/gridpositions1.txt")
data = np.genfromtxt(argv[7]+"/boundarypositions.txt",delimiter=", ",skip_header=(nbounds*nplot),max_rows=nbounds)

fig,ax=plt.subplots()

ax.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
ax.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)

ordered = np.zeros((nbounds,3))
i=0
currentelement = 0
previouselement = 0
while i<nbounds:
    neighbour = int(data[currentelement,-2])
    if neighbour == previouselement:
        neighbour = int(data[currentelement,-1])
    ordered[i,:] = data[neighbour,:3]
    i=i+1
    previouselement=currentelement
    currentelement = neighbour

xs = np.append(ordered[:,0],ordered[0,0])
ys = np.append(ordered[:,1],ordered[0,1])
cs = np.append(ordered[:,2],ordered[0,2])

# Use coloured lines as per https://matplotlib.org/3.1.1/gallery/lines_bars_and_markers/multicolored_line.html
# Create a set of line segments so that we can color them individually
# This creates the points as a N x 1 x 2 array so that we can stack points
# together easily to get the segments. The segments array for line collection
# needs to be (numlines) x (points per line) x 2 (for x and y)
points = np.array([xs, ys]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)


# Create a continuous norm to map from data points to colors
norm = plt.Normalize(cs.min(), cs.max())
lc = LineCollection(segments, cmap='viridis', norm=norm)
# Set the values used for colormapping
lc.set_array(cs)
lc.set_linewidth(2)
line = ax.add_collection(lc)
fig.colorbar(line, ax=ax)

# Use a boundary norm instead
#cmap = ListedColormap(['r', 'g', 'b'])
#norm = BoundaryNorm([-1, -0.5, 0.5, 1], cmap.N)
#lc = LineCollection(segments, cmap=cmap, norm=norm)
#lc.set_array(cs)
#lc.set_linewidth(2)
#line = ax.add_collection(lc)
#fig.colorbar(line, ax=ax)

# ax.plot(xs,ys,c="black")
#ax.scatter(data[:,0],data[:,1],c=data[:,2],s=5,cmap='gray_r',vmin=0,vmax=np.max(data[:,2]),edgecolors='none')
ax.set_xlim([xmin,xmax])
ax.set_ylim([xmin,xmax])
ax.set_aspect('equal')
fig.savefig(argv[7]+"/plot{:04d}.png".format(nplot),bbox_inches='tight',dpi=200)

#os.system("convert output/plot{:04d}.png -shave 268x143 output/plot{:04d}.png".format(nplot,nplot))
