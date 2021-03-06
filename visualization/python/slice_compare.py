# This script calculates (SAM average - CAM average) and plots it in the same way that either average is plotted
#
# Run this from the python directory with:
# >> python slice_compare.py <Name of SAM averaged file in /experiments directory, with extension> <Name of CAM averaged file in /experiments directory, 
#	with extension> <(optional) name of output image, without extension>
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
        
        
#Get the radial helpers
sys.path.insert(0, './../../statistics/')
import radial_helpers

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
	print "call as: \n>>python slice_compare.py <dir, relative to calling directory> <name of SAM file> <name of CAM file> <(optional) name output file, without extension>\n"
	sys.exit(1)
	
	
#define input and output directory
input_dir = sys.argv[1] #'./../../experiments/'
output_dir = input_dir

# Define color of empty cells
empty_color_diff='gray'
empty_color_avgs='white'

#get header
if sys.argv[2]!='':
	
	input_file_sam = input_dir + sys.argv[2];
	input_file_cam = input_dir + sys.argv[3];
	
	fh_sam=open(input_file_sam,'r')
	#fh_cam=open(input_file_cam,'r') # currently the script does not check that the sam and cam data have been obtained from the same raw datafiles.
	
	if len(sys.argv)==5:
		outputname=sys.argv[4]+".pdf" #change format if needed
		print "Plots will be saved as " + outputname + ", cam_" + outputname + ", sam_" + outputname + ', and diffradavg_' + outputname
	else:
		outputname="output_compare.pdf"
		print "Plots will be saved as output_compare.pdf, cam_output_compare.pdf, sam_output_compare.pdf and diffradavg_output_compare.pdf"
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
difference_1d=np.subtract(averages_sam,averages_cam)

#prepare data for plotting
difference = difference_1d.reshape(ny,nz)
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
metadata['Title'] = 'Data plotted (2nd column) = difference ' + sys.argv[1]+sys.argv[2] +' [MINUS] '+sys.argv[1] +sys.argv[3] #input_file_sam, input_file_cam
metadata['Author'] = 'Script used to plot this = slice_compare.py'
metadata['Subject']= 'Info on CAM and SAM avgs plotted in cam_'+outputname+' and sam_'+outputname
metadata['Keywords']= 'If data is an average, first_file='+ str(first_file) + ', last_file=' + str(last_file) + ', stride=' + str(stride)
#metadata['Creator'] = 
#metadata['Producer']=
pdffig.close()


##########################################
#Radial average of comparison plot

#Get radial average points
y_values = []
x_vals = []
num_points = 7*ny;
a = 1.0;  #<--------------------------------------------------------------------THIS SHOULD BE GIVEN AS A PARAMETER
L = ny*a;
L_half = 0.5*L;
R = 0.5*(L-a);

for i in range(0,num_points):
	radius = a + i*(R-a)/(1.0*num_points-1.0);
	x_vals.insert(0,L_half-radius)
	x_vals.append(L_half+radius)
	lengths = radial_helpers.find_seg_quant(a,L,ny,nz,radius,'length')
	warray=radial_helpers.find_weights(lengths,radius)
	y_values.append(   np.dot( np.nan_to_num(difference_1d)  ,  warray )   )
y_vals=np.concatenate((y_values[::-1],y_values))


#Plot
plt.figure(2)
plt.plot(x_vals,y_vals, 'ro', markersize=3)
fig=plt.gcf()
ax=fig.gca()

y_m = y_vals.min()
y_M = y_vals.max()

ax.set_ylim([0.9*y_m,1.1*y_M])
plt.vlines(range(0,ny,1),0.9*y_m,1.1*y_M)
plt.vlines([0.5,ny-0.5],0.9*y_m,1.1*y_M,'g')
ax.set_xlim(ax.get_xlim()[::-1]) 


ax.set_xlabel('Y (cells)',fontweight="bold")
ax.set_ylabel('v_x(SAM)-v_x(CAM)',fontweight="bold") #Latex?		

#save averages plot to disk, prepare metadata object first
pdffig = PdfPages( output_dir + 'diffradavg_' + outputname )
plt.savefig( pdffig, format="pdf")

#add metadata to figure
metadata = pdffig.infodict()
metadata['Title'] = 'Data plotted is radial average of SAM - CAM velocity averages'
metadata['Author'] = 'Script used to plot this = slice_compare.py.py and radial_helpers.py'
metadata['Subject']= ' Datafiles: '+input_file_sam + ' , '+input_file_cam
#metadata['Keywords']= ''
#metadata['Creator'] = 
#metadata['Producer']=

#Save figure
pdffig.close()

#========================================
#auxiliary plot 1/2 (SAM data only)
plt.figure(3)

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
metadata['Title'] = 'Data plotted (2nd column) =' + sys.argv[1]+sys.argv[2] #input_file_sam
metadata['Author'] = 'Script used to plot this = slice_compare.py'
metadata['Subject']= 'Info on corresponding CAM avg plotted in cam_'+outputname+', difference in '+ outputname
metadata['Keywords']= 'If data is an average, first_file='+ str(first_file) + ', last_file=' + str(last_file) + ', stride=' + str(stride)
#metadata['Creator'] = 
#metadata['Producer']=
pdffig_sam.close()


#========================================
#auxiliary plot 2/2 (CAM data only)
plt.figure(4)

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
metadata['Title'] = 'Data plotted (2nd column) =' + sys.argv[1]+sys.argv[3] #input_file_cam
metadata['Author'] = 'Script used to plot this = slice_compare.py'
metadata['Subject']= 'Info on corresponding SAM avg plotted in sam_'+outputname+', difference in '+ outputname
metadata['Keywords']= 'If data is an average, first_file='+ str(first_file) + ', last_file=' + str(last_file) + ', stride=' + str(stride)
#metadata['Creator'] = 
#metadata['Producer']=
pdffig_cam.close()

