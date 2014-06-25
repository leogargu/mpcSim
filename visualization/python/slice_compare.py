# This script calculates SAM average - CAM average and plots it in the same way that eaither average is plotted
#
# Run this from the /python directory with:
# >> python slice_compare.py <Name of SAM averaged file in DATA directory, with extension> <Name of SAM averaged file in DATA directory, with extension> <(optional) name of output image, without extension>
# If the name for the output plot is omitted, the result is saved as "output_compare.png"
#

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sys

#get header
if sys.argv[1]!='':
	
	input_file_sam='./../../DATA/'+sys.argv[1];
	input_file_cam='./../../DATA/'+sys.argv[2];
	
	fh_sam=open(input_file_sam,'r')
	#fh_cam=open(input_file_cam,'r')
	
	if len(sys.argv)==4:
		outputname=sys.argv[2]+".png"
	else:
		outputname="output_compare.png"
		print "Saving plot in output_compare.png\n"
else:
	print "Call as: python slice_compare.py <name of SAM av file> <name of CAM av file> <name of output plot, without the extension>"
	sys.exit(1)
	
header=fh_sam.readline()
fh_sam.close()

#sanity check: compare both headers/assert

# process dimensions of simulation box 
dimensions=header.split('\t')[:3]
nx=int(dimensions[0])
ny=int(dimensions[1])
nz=int(dimensions[2])

#get data for the slice to be plotted, save in pandas dataframe
data_sam = pd.read_table(input_file_sam,sep='\t',skiprows=1,skipinitialspace=True,header=None)
data_cam = pd.read_table(input_file_cam,sep='\t',skiprows=1,skipinitialspace=True,header=None)

#separate columns
averages_sam=data_sam.values[:,1]
averages_cam=data_cam.values[:,0]

#calculate difference
difference=np.subtract(averages_sam,averages_cam)

print difference


#prepare data for plotting
difference=difference.reshape(ny,nz)



#plot difference
im_difference=plt.imshow(difference,interpolation='bicubic')#or bilinear, nearest, bicubic
plt.colorbar(im_averages,orientation='vertical')

#save averages plot to disk
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