# This script calculates (SAM average - CAM average) and plots it in the same way that either average is plotted
#
# Run this from the python directory with:
# >> python slice_compare.py <Name of SAM averaged file in /experiments directory, with extension> <Name of CAM averaged file in /experiments directory, with extension> <(optional) name of output image, without extension>
# If the name for the output plot is omitted, the result is saved as "output_compare.pdf"
#

import pandas as pd
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
import matplotlib.pyplot as plt
from matplotlib.colors import from_levels_and_colors
from matplotlib import ticker
from matplotlib.backends.backend_pdf import PdfPages
#On PdfPages:
#http://stream.princeton.edu/AWCM/LIBRARIES/matplotlib-1.3.0/lib/matplotlib/backends/backend_pdf.py        
        
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
               
# check number of calling arguments
if len(sys.argv)<3:
	print "call as: \n>>python slice_compare.py <name of SAM file in /experiments> <name of CAM file in /experiments> <(optional) name output file, without extension>\n"
	sys.exit(1)
	
	
#define input and output directory
input_dir = './../../experiments/'
output_dir = input_dir

# Define color of empty cells
empty_color_diff='gray'
empty_color_avgs='white'

#get header
if sys.argv[1]!='':
	
	input_file_sam = input_dir + sys.argv[1];
	input_file_cam = input_dir + sys.argv[2];
	
	fh_sam=open(input_file_sam,'r')
	#fh_cam=open(input_file_cam,'r') # currently the script does not check that the sam and cam data have been obtained from the same raw datafiles.
	
	if len(sys.argv)==4:
		outputname=sys.argv[3]+".pdf" #change format if needed
		print "Plots will be saved as " + outputname + ", cam_" + outputname + ", and sam_" + outputname
	else:
		outputname="output_compare.pdf"
		print "Plots will be saved as output_compare.pdf, cam_output_compare.pdf and sam_output_compare.pdf"
	print "In ", output_dir
else:
	print "Call as: python slice_compare.py <name of SAM av file> <name of CAM av file> <name of output plot, without the extension>"
	sys.exit(1)
	
header=fh_sam.readline()
fh_sam.close()

#sanity check: compare both headers/assert

# process dimensions of simulation box 
dimensions=header.split()
nx=int(dimensions[0])
ny=int(dimensions[1])
nz=int(dimensions[2])
idx_first_cell = int(dimensions[3])	#global index of first cell included in datafile
idx_last_cell = int(dimensions[4])	#global index of last cell included in datafile
first_file = int(dimensions[5])
stride = int(dimensions[6])
last_file = int(dimensions[7])


#get data for the slice to be plotted, save in pandas dataframe
data_sam = pd.read_table(input_file_sam,sep='\t',skiprows=1,skipinitialspace=True,header=None)
data_cam = pd.read_table(input_file_cam,sep='\t',skiprows=1,skipinitialspace=True,header=None) #if header is skipped, there is only one column that starts with an empty row: read_table doesn't know how to parse this.

#separate columns
averages_sam = data_sam.values[:,1]
averages_cam = data_cam.values[:,1]

#find min and max over sam and cam combined to adjust auxiliary CAM/SAM plots
my_vmin=np.nanmin(averages_sam)
if np.nanmin(averages_cam) < my_vmin:
	my_vmin = np.nanmin(averages_cam)
my_vmax=np.nanmax(averages_sam)
if np.nanmax(averages_cam) > my_vmax:
	my_vmax = np.nanmax(averages_cam)	

#calculate difference
difference=np.subtract(averages_sam,averages_cam)

#prepare data for plotting
difference = difference.reshape(ny,nz)
averages_cam = averages_cam.reshape(ny,nz)
averages_sam = averages_sam.reshape(ny,nz)

##########################################
#Plot comparison plot
num_levels = 50
vmin, midpoint, vmax = np.nanmin(difference), 0.0,  np.nanmax(difference)
levels = np.linspace(vmin, vmax, num_levels)
midp = np.mean(np.c_[levels[:-1], levels[1:]], axis=1)

vals = np.interp(midp, [vmin, midpoint, vmax], [0, 0.5, 1])
colors = plt.cm.seismic(vals)
cmap, norm = from_levels_and_colors(levels, colors)  #won't use the norm

masked_array = np.ma.array(difference, mask=np.isnan(difference))
cmap.set_bad(color='gray')

fig, ax = plt.subplots()
im = ax.imshow(masked_array, vmin=vmin, vmax=vmax,cmap=cmap,interpolation=interpolation_setting,extent=[0,nz,0,ny])
cb=fig.colorbar(im,orientation='vertical')
cb.locator = ticker.MaxNLocator(nbins=num_levels/3,symmetric=True,prune=None)
cb.update_ticks()

plot_refs(plt,ny,nz)


#save averages plot to disk
pdffig = PdfPages( output_dir + outputname )
plt.savefig( pdffig, format="pdf")

#add metadata to figure
metadata = pdffig.infodict()
metadata['Title'] = 'Data plotted (2nd column) = difference ' + sys.argv[1] +' [MINUS] '+ sys.argv[2] #input_file_sam, input_file_cam
metadata['Author'] = 'Script used to plot this = slice_compare.py'
metadata['Subject']= 'Info on CAM and SAM avgs plotted in cam_'+outputname+' and sam_'+outputname
metadata['Keywords']= 'If data is an average, first_file='+ str(first_file) + ', last_file=' + str(last_file) + ', stride=' + str(stride)
#metadata['Creator'] = 
#metadata['Producer']=
pdffig.close()




########################remove this

original_cmap = plt.get_cmap('seismic')
remapped_cmap = col.remappedColorMap(original_cmap, data=difference)

#Plots:
#(1) Original Cmap
plt.subplot(2, 1, 1)
test_plot2=plt.imshow(difference,cmap=plt.get_cmap('seismic') ,interpolation='bicubic')
plt.title("Original Cmap")
divider = make_axes_locatable(plt.gca())
plt.colorbar(test_plot2,orientation='vertical', cax=divider.append_axes("right", size="5%", pad=0.05))

#(2) Remapped Cmap
plt.subplot(2, 1, 2)
test_plot1=plt.imshow(difference,cmap=remapped_cmap ,interpolation='nearest')
plt.title("Remapped Cmap")
divider = make_axes_locatable(plt.gca())
plt.colorbar(test_plot1,orientation='vertical', cax=divider.append_axes("right", size="5%", pad=0.05))

plt.subplots_adjust(hspace=.3)

#Saving plot
plt.savefig("test_fig.png",bbox_inches='tight')


##########################################
#auxiliary plot 1/2 (SAM data only)
plt.figure(2)

# plot empty cells in a pre-determined color
my_cmap=plt.get_cmap('jet') #Also interesting: jet, spectral, rainbow 
masked_array = np.ma.array(averages_sam, mask=np.isnan(averages_sam))
my_cmap.set_bad(color=empty_color_avgs)

im_sam=plt.imshow(averages_sam,cmap=my_cmap,vmin=my_vmin,vmax=my_vmax,interpolation=interpolation_setting,extent=[0,nz,0,ny])#or bilinear, nearest, bicubic; plt.get_cmap('seismic')
plt.colorbar(im_sam,orientation='vertical')

plot_refs(plt,ny,nz)

pdffig_sam = PdfPages( output_dir + "sam_" + outputname )

plt.savefig( pdffig_sam, format="pdf" )

#add metadata to figure
metadata = pdffig_sam.infodict()
metadata['Title'] = 'Data plotted (2nd column) =' + sys.argv[1] #input_file_sam
metadata['Author'] = 'Script used to plot this = slice_compare.py'
metadata['Subject']= 'Info on corresponding CAM avg plotted in cam_'+outputname+', difference in '+ outputname
metadata['Keywords']= 'If data is an average, first_file='+ str(first_file) + ', last_file=' + str(last_file) + ', stride=' + str(stride)
#metadata['Creator'] = 
#metadata['Producer']=
pdffig_sam.close()


#========================================
#auxiliary plot 2/2 (CAM data only)
plt.figure(3)

# plot empty cells in a pre-determined color
my_cmap=plt.get_cmap('jet') #Also interesting: jet, spectral, rainbow 
masked_array = np.ma.array(averages_cam, mask=np.isnan(averages_cam))
my_cmap.set_bad(color=empty_color_avgs)

im_cam=plt.imshow(averages_cam,cmap=my_cmap,vmin=my_vmin,vmax=my_vmax,interpolation=interpolation_setting,extent=[0,nz,0,ny])#or bilinear, nearest, bicubic; plt.get_cmap('seismic')
plt.colorbar(im_cam,orientation='vertical')

plot_refs(plt,ny,nz)


pdffig_cam = PdfPages( output_dir + "cam_" + outputname )

plt.savefig( pdffig_cam , format="pdf")

#add metadata to figure
metadata = pdffig_cam.infodict()
metadata['Title'] = 'Data plotted (2nd column) =' + sys.argv[2] #input_file_cam
metadata['Author'] = 'Script used to plot this = slice_compare.py'
metadata['Subject']= 'Info on corresponding SAM avg plotted in sam_'+outputname+', difference in '+ outputname
metadata['Keywords']= 'If data is an average, first_file='+ str(first_file) + ', last_file=' + str(last_file) + ', stride=' + str(stride)
#metadata['Creator'] = 
#metadata['Producer']=
pdffig_cam.close()

