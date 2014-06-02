# Before running this script, calculate velprof_CAM_averaged.dat by running this from /statistics:
# >> ./average CAM velprof <first file index> <last file index> 1.0 <verbose>
#
# Then, run this from the /python directory with:
# >> python velprofile.py "Name of file in DATA directory"
#
# This generates an image, velprofile.png, of the cross-sectional velocity profile, for example.


import matplotlib.pyplot as plt
import matplotlib.image as mpimg
#from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sys



#get data for the slice to be plotted
data=np.genfromtxt('./../../DATA/equilibration.dat',delimiter='\t')

#plot momentum
plt.plot(data)
plt.legend(["Jx","Jy","Jz"],loc=7)
plt.xlabel("timestep")
plt.ylabel("Total linear momentum J")
plt.title("Equilibration phase")


# Save it to disk
plt.savefig("./momentum.png") #how to save eps?

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
