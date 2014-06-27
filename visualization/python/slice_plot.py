# Before running this script, calculate velprof_SAM_averaged.dat or SAMconverted_velprof_SAM_averaged.dat by running this from /statistics:
# >> ./average SAM/CAM <name of SAM/CAM-formatted file> <first file index> <last file index> <factor> <verbose>
#
#
# Then, run this from the /python directory with:
# >> python slice_plot.py <Name of file in /experiments directory, with extension> <(optional) Name of output plot (without the extension)>
# If the name for the output plot is omitted, the result is saved as "output.png"
# a plot containing information on the number of samples used (for SAM averages) or number of particles across samples (for CAM averages)
# is also plotted and saved in "output_samples.png" , or the corresponding given name with the suffix "_samples.png".
#
# This script generates an image, velprofile.png, of the cross-sectional velocity profile, for example.

import pandas as pd
import matplotlib.pyplot as plt
#import matplotlib.image as mpimg
#from mpl_toolkits.mplot3d import Axes3D
import pylab as pl
import numpy as np
import sys
from matplotlib.backends.backend_pdf import PdfPages
#On PdfPages:
#http://stream.princeton.edu/AWCM/LIBRARIES/matplotlib-1.3.0/lib/matplotlib/backends/backend_pdf.py

#define input and output directory
input_dir = './../../experiments/'
output_dir = input_dir

#check health of input arguments and get header
if len(sys.argv)>=2 and sys.argv[1]!='':
	
	input_file=input_dir + sys.argv[1];
	
	fh=open(input_file,'r')
	if len(sys.argv)==3:
		outputname=sys.argv[2]+".pdf" #change the extension if required
		outputname_samples=sys.argv[2]+"_samples.pdf"
	else:
		outputname="output.pdf"
		outputname_samples="output_samples.pdf"
		print "Plots will be saved in output.pdf and output_samples.pdf\n"
else:
	print "Call as: \n>>python slice_plot.py <name of av file to plot, in /experiments> <(optional) name of output plot, without the extension>"
	sys.exit(1)
	
header=fh.readline()
fh.close()

# process dimensions of simulation box
dimensions=header.split()
nx=int(dimensions[0])
ny=int(dimensions[1])
nz=int(dimensions[2])
idx_first_cell = int(dimensions[3])
idx_last_cell = int(dimensions[4])
first_file = int(dimensions[5])
last_file = int(dimensions[6])

#get data for the slice to be plotted, save in pandas dataframe
data = pd.read_table(input_file,sep='\t',skiprows=1,skipinitialspace=True,header=None)

#separate columns
averages=data.values[:,1]
samples=data.values[:,0]

data_min = np.nanmin(averages)
data_max = np.nanmax(averages)

#in here, find the max and min of samples, and print to screen
print "Maximum number of samples(SAM)/particles(CAM) used: ", int(max(samples))
samples_aux=sorted(set(samples)) #remove duplicated entries, return a sorted list
print "Minimum number of samples(SAM)/particles(CAM) used: ", int(samples_aux[1])

#prepare data for plotting
samples=samples.reshape(ny,nz)
averages=averages.reshape(ny,nz)

##########################################

#plot samples statistics
samples[samples==0] = np.nan

im_samples=plt.imshow(samples,interpolation='nearest')
plt.colorbar(im_samples, orientation='vertical')

# Save sample statistics plot to disk
pdffig_samples = PdfPages( output_dir + outputname_samples )
#plt.savefig( output_dir + outputname_samples )
plt.savefig( pdffig_samples, format="pdf" )

#add metadata to figure
metadata = pdffig_samples.infodict()
metadata['Title'] = 'Data plotted =' + input_file 
metadata['Author'] = 'Script used to plot this = slice_plot.py'
metadata['Subject']= 'Info on averages plotted in '+outputname
metadata['Keywords']= 'If data is an average, first_file='+ str(first_file) + ', last_file=' + str(last_file)
#metadata['Creator'] = 
#metadata['Producer']=
pdffig_samples.close()

##########################################

#plot averages
plt.figure(2)
im_averages=plt.imshow(averages, vmin=data_min,vmax=data_max,interpolation='bicubic')#or bilinear, nearest, bicubic
plt.colorbar(im_averages,orientation='vertical')

#save averages plot to disk, prepare metadata object first
pdffig = PdfPages( output_dir + outputname )
plt.savefig( pdffig, format="pdf")#output_dir + outputname )

#add metadata to figure
metadata = pdffig.infodict()
metadata['Title'] = 'Data plotted =' + input_file 
metadata['Author'] = 'Script used to plot this = slice_plot.py'
metadata['Subject']= 'Info on samples(SAM)/particles(CAM) plotted in '+outputname_samples
metadata['Keywords']= 'If data is an average, first_file='+ str(first_file) + ', last_file=' + str(last_file)
#metadata['Creator'] = 
#metadata['Producer']=

#Note that the time is 1 hour off in the metadata of the pdf if opened with okular, even though it is alright in the metadata dict
#Evince displays the metadata date and time correctly

#Metadata fields in pdf documents are restricted to set fields. Those that accept text strings are: Title, Author, Subject,Keywords,Creator,Producer. 
#In addition, there is CReationDate and ModDate

pdffig.close()

##########################################

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