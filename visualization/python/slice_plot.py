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
import numpy as np
#from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
import matplotlib.pyplot as plt
import matplotlib

from matplotlib.colors import from_levels_and_colors
from matplotlib import ticker
#import pylab as pl
from matplotlib.backends.backend_pdf import PdfPages
#On PdfPages:
#http://stream.princeton.edu/AWCM/LIBRARIES/matplotlib-1.3.0/lib/matplotlib/backends/backend_pdf.py
#from slice_compare import plot_refs (something didn't work well: error message is taken from slice_compare instead of slice_plot)

# Auxiliary functions
def plot_refs(plt,ny,nz):
	#Beautify the axes and ticks, get them to conform to standard reference
	ax=plt.gca()
	ax.arrow(0, ny, nz, 0, fc='k', ec='k',lw = 2,head_width=0.1, head_length=0.2, length_includes_head= True, clip_on = False)
	ax.arrow(0, nz, 0., -nz, fc='k', ec='k',lw = 2, head_width=0.1, head_length=0.2, length_includes_head= True, clip_on = False)
	ax.get_xaxis().set_ticks([])
	ax.get_yaxis().set_ticks([])
	ax.set_ylabel('Y',fontweight="bold")
	ax.xaxis.set_label_position('top') 
	ax.set_xlabel('Z',fontweight="bold")
	
	# Plot boundary circunference
	boundary = plt.Circle((ny*0.5,nz*0.5),radius=(ny-1)*0.5,linewidth=1,color='k',fill=False)
	fig = plt.gcf()
	fig.gca().add_artist(boundary)	
	return

# Interpolation setting
interpolation_setting='nearest' #also possible: bilinear, bicubic
 

#define input and output directory
input_dir = './../../experiments/'
output_dir = input_dir

# Define color of empty cells
empty_color_stats='black'
empty_color_avgs='white'

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

is_snapshot = False
is_average = False

# process dimensions of simulation box
dimensions=header.split()
if len(dimensions) == 6:
	print "Plotting a snapshot..."
	is_snapshot = True
elif len(dimensions) == 8:
	print "Plotting an average..."
	is_average = True
else:
	print "WARNING: Check datafile is not corrupt (not plotting snapshot nor average)"
	sys.exit(1)
	


nx=int(dimensions[0])
ny=int(dimensions[1])
nz=int(dimensions[2])
idx_first_cell = int(dimensions[3])	#global index of the first cell included in the datafile
idx_last_cell = int(dimensions[4])	#global index of the last cell included in the datafile
#if the datafile is an average, header is: nx ny nz first_cell last_cell first_file stride last_file
if is_average:
	first_file = int(dimensions[5])		#first datafile number used in the average that generated the datafile plotted here 
	stride = int(dimensions[6])		#stride used when calculating the average datafile plotted here
	last_file = int(dimensions[7])		#last datafile number used in the average that generated the datafile plotted here
#if the datafile is a snapshot, header is: nx ny nz first_cell last_cell timestep
if is_snapshot:
	timestep = int(dimensions[5])

#get data for the slice to be plotted, save in pandas dataframe
data = pd.read_table(input_file,sep='\t',skiprows=1,skipinitialspace=True,header=None)

#separate columns
averages=data.values[:,1]
samples=data.values[:,0]

data_min = np.nanmin(averages)
data_max = np.nanmax(averages)

#in here, find the max and min of samples, and print to screen
if is_average:
	print "Maximum number of samples(SAM)/particles(CAM) used: ", int(max(samples))
	samples_aux=sorted(set(samples)) #remove duplicated entries, return a sorted list
	print "Minimum number of samples(SAM)/particles(CAM) used: ", int(samples_aux[1])

#prepare data for plotting
samples=samples.reshape(ny,nz)
averages=averages.reshape(ny,nz)

##########################################
#The following two sections could be rewritten as 1-2 functions?
##########################################
#plot samples statistics
samples[samples==0] = np.nan

my_cmap=plt.get_cmap('Greens') #Also interesting: jet, spectral, rainbow 
masked_array = np.ma.array(samples, mask=np.isnan(samples))
my_cmap.set_bad(color=empty_color_stats)

im_samples=plt.imshow(samples,cmap=my_cmap,interpolation='nearest',extent=[0,nz,0,ny])
plt.colorbar(im_samples, orientation='vertical')


plot_refs(plt,ny,nz)

# Save sample statistics plot to disk
pdffig_samples = PdfPages( output_dir + outputname_samples )
plt.savefig( pdffig_samples, format="pdf" )

#add metadata to figure
metadata = pdffig_samples.infodict()
metadata['Title'] = 'Data plotted =' + sys.argv[1] #input_file 
metadata['Author'] = 'Script used to plot this = slice_plot.py'
metadata['Subject']= 'Info on averages plotted in '+outputname
if is_average:
	metadata['Keywords']= 'If data is an average, first_file='+ str(first_file) + ', last_file=' + str(last_file) + ', stride=' + str(stride)
if is_snapshot:
	metadata['Keywords']='Data is a snapshot taken at timestep='+str(timestep)
#metadata['Creator'] = 
#metadata['Producer']=
pdffig_samples.close()


##########################################
#plot averages or snapshot
plt.figure(2)
my_cmap=plt.get_cmap('jet') #Also interesting: jet, spectral, rainbow 
masked_array = np.ma.array(averages, mask=np.isnan(averages))
my_cmap.set_bad(color=empty_color_avgs)

im_averages=plt.imshow(averages,cmap=my_cmap,vmin=data_min,vmax=data_max,interpolation=interpolation_setting, extent=[0,nz,0,ny])#or bilinear, nearest, bicubic
plt.colorbar(im_averages,orientation='vertical')

#Beautify the axes and ticks, get them to conform to standard reference
plot_refs(plt,ny,nz)


#save averages plot to disk, prepare metadata object first
pdffig = PdfPages( output_dir + outputname )
plt.savefig( pdffig, format="pdf")#output_dir + outputname )

#add metadata to figure
metadata = pdffig.infodict()
metadata['Title'] = 'Data plotted =' + sys.argv[1] #input_file 
metadata['Author'] = 'Script used to plot this = slice_plot.py'
metadata['Subject']= 'Info on samples(SAM)/particles(CAM) plotted in '+outputname_samples
if is_average:
	metadata['Keywords']= 'Data is an average, first_file='+ str(first_file) + ', last_file=' + str(last_file) + ', stride=' + str(stride)
if is_snapshot:
	metadata['Keywords']= 'Data is a snapshot taken at timestep='+str(timestep)
#metadata['Creator'] = 
#metadata['Producer']=

#Note that the time is 1 hour off in the metadata of the pdf if opened with okular, even though it is alright in the metadata dict
#Evince displays the metadata date and time correctly

#Metadata fields in pdf documents are restricted to set fields. Those that accept text strings are: Title, Author, Subject,Keywords,Creator,Producer. 
#In addition, there is CReationDate and ModDate

pdffig.close()