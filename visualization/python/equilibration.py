# Run this from the /python directory with:
# >> python equilibration.py 
#
# This generates an image, momentum.png, of the evolution of total mmentum at each timestep in the simulation


#NOT TESTED after changing directories

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
#from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sys

#define input and output directory
input_dir = './../../experiments/'
output_dir = input_dir

#get data for the slice to be plotted
data=np.genfromtxt(input_dir+'equilibration.dat',delimiter='\t')

#plot momentum
plt.plot(data)
plt.legend(["Jx","Jy","Jz"],loc=7)
plt.xlabel("timestep")
plt.ylabel("Total linear momentum J")
plt.title("Equilibration phase")


# Save it to disk
plt.savefig(output_dir+"momentum.png") #how to save eps?

#add functionality to plot kinetic energy of the system  here

#--------------------
#plot 3d wireframe
#--------------------
#X = np.arange(0, nz)
#Y = np.arange(0, ny)
#X, Y = np.meshgrid(X, Y)
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#surf = ax.plot_wireframe(Y, X, data)
#ax.view_init(30, 0)
