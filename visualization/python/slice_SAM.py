# Before running this script, calculate velprof_SAM_averaged.dat or SAMconverted_velprof_SAM_averaged.dat by running this from /statistics:
# >> ./average SAM <name of SAM-formatted file> <first file index> <last file index> 1.0 <verbose>
#
#
# Then, run this from the /python directory with:
# >> python slice_SAM.py <Name of file in DATA directory, with extension> <Name of output plot (without the extension)>
# If the name for the output plot is omitted, the result is saved as "output.png"
# a plot containing information on the number of samples used is also plotted and saved in "output_samples.png", or the corresponding given name 
# with the suffix "_samples"
#
# This generates an image, velprofile.png, of the cross-sectional velocity profile, for example.

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sys

#get header
if sys.argv[1]!='':
	
	input_file='./../../DATA/'+sys.argv[1];
	
	fh=open(input_file,'r')
	if len(sys.argv)==3:
		outputname=sys.argv[2]+".png"
		outputname_samples=sys.argv[2]+"_samples.png"
	else:
		outputname="output.png"
		outputname_samples="output_samples.png"
		print "Saving plot in output.png and output_samples.png\n"
else:
	print "Call as: python slice_SAM.py <name of av file to plot> <name of output plot, without the extension>"
	sys.exit(1)
	
header=fh.readline()
fh.close()

# process dimensions of simulation box
dimensions=header.split('\t')[:3]
nx=int(dimensions[0])
ny=int(dimensions[1])
nz=int(dimensions[2])

#get data for the slice to be plotted, save in pandas dataframe
data = pd.read_table(input_file,sep='\t',skiprows=1,skipinitialspace=True,header=None)

#separate columns
averages=data.values[:,1]
samples=data.values[:,0]


#in here, find the max and min of samples, and print to screen
print "Maximum number of samples used: ", int(max(samples))
samples_aux=sorted(set(samples)) #remove duplicated entries, return a sorted list
print "Minimum number of samples used: ", int(samples_aux[1])

#prepare data for plotting
samples=samples.reshape(ny,nz)
averages=averages.reshape(ny,nz)

#plot samples statistics
im_samples=plt.imshow(samples,interpolation='nearest')
plt.colorbar(im_samples, orientation='vertical')

# Save sample statistics plot to disk
plt.savefig("./"+outputname_samples)#how to save eps?


#plot averages
plt.figure(2)
im_averages=plt.imshow(averages, vmin=0,vmax=6,interpolation='nearest')#or bilinear, nearest, bicubic
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
