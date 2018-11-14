import numpy as np
import matplotlib.pyplot as plt
import os

data = np.genfromtxt("output/boundarypositions.txt",delimiter=", ")
nbounds = np.genfromtxt("output/nbounds.txt",delimiter=", ",dtype=int)
#nbounds.astype(int)
Nc = 3
xmax=10
npoints = nbounds[-1]
drawn = 0
for step in range(nbounds.shape[0]):
    print("{:02d}/{:02d}".format((step+1),nbounds.shape[0]))
    fig, ax = plt.subplots(figsize=(4,4))
    print("{}:{}".format(drawn,drawn+nbounds[step]))

    ax.scatter(data[drawn:drawn+nbounds[step],0],data[drawn:drawn+nbounds[step],1])
    drawn = drawn+nbounds[step]
    ax.set_xlim([-xmax/10,xmax/10])
    ax.set_ylim([-xmax/10,xmax/10])
    ax.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
    ax.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)
    fig.savefig("output/test{:02d}".format(step))
    plt.close()
os.system("convert -delay 10 -loop 0 output/*.png output/animated.gif;")
# rm output/*.png")
