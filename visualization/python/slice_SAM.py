# Before running this script, calculate velprof_SAM_averaged.dat by running this from /statistics:
# >> ./average SAM velprof <first file index> <last file index> 1.0 <verbose>
#
# Then, run this from the /python directory with:
# >> python slice_SAM.py <Name of file in DATA directory, with extension> <Name of output plot (without the extension)>
# If the name for the output plot is omitted, the result is saved as "output.png"
#
# This generates an image, velprofile.png, of the cross-sectional velocity profile, for example.

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.mplot3d import Axes3D
#import numpy as np
import sys

#get header
if sys.argv[1]!='':
	if sys.argv[1]=='test':
		input_file='./tests/velprof_SAM_test2.dat';
	else:		
		input_file='./../../DATA/'+sys.argv[1];
	
	fh=open(input_file,'r')
	if len(sys.argv)==3:
		outputname=sys.argv[2]+".png"
	else:
		outputname="output.png"
		print "Saving plot in output.png\n"
else:
	print "Call as: python slice_SAM.py <name of av file to plot> <name of output plot,without the extension>"
	sys.exit(1)
	
header=fh.readline()
fh.close()

# process dimensions of simulation box
dimensions=header.split('\t')[:3]
nx=int(dimensions[0])
ny=int(dimensions[1])
nz=int(dimensions[2])

#get data for the slice to be plotted
#data=np.genfromtxt(input_file,skip_header=1,delimiter='\t',usecols=(0,1))

data = pd.read_table('data',sep='\t',lineterminator='\n',skiprows=1,header=None)
print data

#averages=data[:,1]
#samples=data[:,0]

#print averages
#print samples

#in here, find the max and min of samples, and print to screen
#print "Maximum number of samples used: ", int(max(samples))," <=total number of samples used in the average"
#samples_aux=sorted(set(samples)) #remove duplicated entries, return a sorted list
#print "Minimum number of samples used: ", int(samples_aux[1])

#prepare data for plotting
#samples=samples.reshape(ny,nz)
#averages=averages.reshape(ny,nz)


#plot
#im=plt.imshow(averages,interpolation='nearest')#or bilinear, nearest, bicubic
#plt.colorbar(im, orientation='vertical')

# Save it to disk
#plt.savefig("./"+outputname)#how to save eps?


#Show the number of samples used
#im_samples=plt.imshow(samples,interpolation='nearest')#or bilinear, nearest, bicubic
#plt.colorbar(im_samples, orientation='vertical')
#change above: two vertical bars appear. what is the handle for the image?

# Save it to disk
#plt.savefig("./samples_"+outputname)#how to save eps?



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
