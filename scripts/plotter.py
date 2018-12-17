import numpy as np
import matplotlib.pyplot as plt
import os

os.system("rm output/*.png;rm output/animatedboundary.gif;")
data = np.genfromtxt("output/boundarypositions.txt",delimiter=", ")
nbounds = np.genfromtxt("output/nbounds.txt",delimiter=", ",dtype=int)
xmax=0.5
drawn = 0
fig, ax = plt.subplots(figsize=(4,4))
for step in range(nbounds.shape[0]):
    if step%10 < 1:
        print("{:02d}/{:02d}".format((step+1),nbounds.shape[0]))
        ax.set_xlim([-xmax,xmax])
        ax.set_ylim([-xmax,xmax])
        ax.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
        ax.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)
        ax.plot(np.append(data[drawn:drawn+nbounds[step],0],data[drawn,0]),np.append(data[drawn:drawn+nbounds[step],1],data[drawn,1]))
        #ax.scatter(data[drawn:drawn+nbounds[step],0],data[drawn:drawn+nbounds[step],1])
        fig.savefig("output/test{:05d}".format(step),bbox_inches='tight',padding_inches=0,dpi=200)
        ax.cla()
    drawn = drawn+nbounds[step]
    plt.close()
os.system("convert -delay 10 -loop 0 output/*.png output/animatedboundary.gif;")
# rm output/*.png")
