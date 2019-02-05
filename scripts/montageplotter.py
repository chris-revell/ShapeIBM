
import import matplotlib.pyplot as plt
import numpy as np
data = np.genfromtxt("output/volume.txt")

fig,ax = plt.subplots()
for i in range(np.shape(data)[0]):
    ax.set_ylim([0,np.max(data[:,1])])
    ax.set_xlim([0,np.max(data[:,0])])
    ax.set_ylabel("Volume /$\mu m^2$")
    ax.set_xlabel("Time /s")
    ax.plot(data[:i,0],data[:i,1])
    fig.savefig("output/volume{:04d}.png".format(i),bbox_inches="tight",padding_inches=0)
    plt.cla()

for i in range(np.shape(data)[0]):
    os.system("montage output/velocitytest{:04d}.png output/volume{:04d}.png -tile 2x1 -geometry 900x output/montage{:04d}.png".format(i,i,i))



os.system("convert -delay 10 -loop 0 output/montage*.png output/montageanimated.gif")
