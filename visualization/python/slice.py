# Before running this script, calculate velprof_CAM_averaged.dat by running this from /statistics:
# >> ./average CAM velprof <first file index> <last file index> 1.0 <verbose>
#
# Then, run this from the /python directory with:
# >> python velprofile.py <Name of file in DATA directory> <Name of output plot (without the extension)>
# If the name for the output plot is omitted, the result is saved as "output.png"
#
# This generates an image, velprofile.png, of the cross-sectional velocity profile, for example.


import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sys

#get header
if sys.argv[1]!='':
	fh=open('./../../DATA/'+sys.argv[1],'r')
	if len(sys.argv)==3:
		outputname=sys.argv[2]+".png"
	else:
		outputname="output.png"
		print "Saving plot in output.png\n"
else:
	print "Call as: python test.py <name of av file to plot> <name of output plot,without the extension>"
	sys.exit(1)
	
header=fh.readline()
fh.close()

# process dimensions of simulation box
dimensions=header.split('\t')[:3]
nx=int(dimensions[0])
ny=int(dimensions[1])
nz=int(dimensions[2])


#get data for the slice to be plotted
data=np.genfromtxt('./../../DATA/'+sys.argv[1],skip_header=1,delimiter='\n')


#prepare data for plotting
data=data.reshape(ny,nz)


#plot
im=plt.imshow(data,interpolation='nearest')#or bilinear, nearest, bicubic
plt.colorbar(im, orientation='vertical')

# Save it to disk
plt.savefig("./"+outputname)#how to save eps?

#how to add axis labels?

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