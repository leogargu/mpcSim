# This script calculates SAM average - CAM average and plots it in the same way that eaither average is plotted
#
# Run this from the /python directory with:
# >> python slice_compare.py <Name of SAM averaged file in DATA directory, with extension> <Name of SAM averaged file in DATA directory, with extension> <(optional) name of output image, without extension>
# If the name for the output plot is omitted, the result is saved as "output_compare.png"
#

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sys

################################################################################################
#This function is by Paul H, modified by TheChymera and comes from here:
#http://stackoverflow.com/questions/7404116/defining-the-midpoint-of-a-colormap-in-matplotlib
#https://github.com/TheChymera/chr-helpers/blob/d05eec9e42ab8c91ceb4b4dcc9405d38b7aed675/chr_matplotlib.py
################################################################################################
from mpl_toolkits.axes_grid1 import AxesGrid

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    -----
    cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and 0.5; if your dataset mean is negative you should leave 
          this at 0.0, otherwise to (vmax-abs(vmin))/(2*vmax) 
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0; usually the
          optimal value is abs(vmin)/(vmax+abs(vmin)) 
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          0.5 and 1.0; if your dataset mean is positive you should leave 
          this at 1.0, otherwise to (abs(vmin)-vmax)/(2*abs(vmin)) 
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.hstack([
        np.linspace(start, 0.5, 128, endpoint=False), 
        np.linspace(0.5, stop, 129)
    ])

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap
    
################################################################################################
        
    


#get header
if sys.argv[1]!='':
	
	input_file_sam='./../../DATA/'+sys.argv[1];
	input_file_cam='./../../DATA/'+sys.argv[2];
	
	fh_sam=open(input_file_sam,'r')
	#fh_cam=open(input_file_cam,'r')
	
	if len(sys.argv)==4:
		outputname=sys.argv[3]+".png"
		print "I will save plots as " + outputname + ", cam_" + outputname + ", and sam_" + outputname +"\n"

	else:
		outputname="output_compare.png"
		print "I will save plots as output_compare.png, cam_output_compare.png and sam_output_compare.png\n"
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
data_cam = pd.read_table(input_file_cam,skipinitialspace=True,header=None) #if header is skipped, there is only one column that starts with an empty row: read_table doesn't know how to parse this.

#separate columns
data_cam_clean=data_cam.values[:,0]
averages_cam=np.delete(data_cam_clean,0)
averages_sam=data_sam.values[:,1]

my_vmin=np.nanmin(averages_sam)
if np.nanmin(averages_cam) < my_vmin:
	my_vmin = np.nanmin(averages_cam)
my_vmax=np.nanmax(averages_sam)
if np.nanmax(averages_cam) > my_vmax:
	my_vmax = np.nanmax(averages_cam)	


#calculate difference
difference=np.subtract(averages_sam,averages_cam)

#the mean of the dataset is used to offset the color map later
data_mean = np.nanmean(difference)
v_min=np.nanmin(difference)    
v_max=np.nanmax(difference)

#prepare data for plotting
difference = difference.reshape(ny,nz)

averages_cam = averages_cam.reshape(ny,nz)
averages_sam = averages_sam.reshape(ny,nz)

#plot difference
mycmap_midpoint = abs(v_min)/(v_max+abs(v_min))

if data_mean<0:
	mycmap_start=0.0
else:
	mycmap_start=0.5*(v_max-abs(v_min))/v_max
	
if data_mean >0:
	mycmap_stop=1.0
else:
	mycmap_stop=0.5*(abs(v_min)-v_max)/abs(v_min)
	
assert mycmap_start<mycmap_midpoint
assert mycmap_midpoint<mycmap_stop
assert 0.5<=mycmap_stop<=1.0
assert 0.0<=mycmap_start<=0.5

orig_cmap = matplotlib.cm.seismic
shifted_cmap = shiftedColorMap(orig_cmap, midpoint=mycmap_midpoint, name='shifted')
shrunk_cmap = shiftedColorMap(orig_cmap, start=mycmap_start, midpoint=mycmap_midpoint, stop=mycmap_stop, name='shrunk')

im_difference=plt.imshow(difference,cmap=shrunk_cmap,interpolation='bicubic')#or bilinear, nearest, bicubic; plt.get_cmap('seismic')
plt.colorbar(im_difference,orientation='vertical')

#save averages plot to disk
plt.savefig("./"+outputname)#how to save eps?

#auxiliary plots
plt.figure(2)
im_sam=plt.imshow(averages_sam,vmin=my_vmin,vmax=my_vmax,interpolation='bicubic')#or bilinear, nearest, bicubic; plt.get_cmap('seismic')
plt.colorbar(im_sam,orientation='vertical')
plt.savefig("./sam_"+outputname)#how to save eps?

plt.figure(3)
im_cam=plt.imshow(averages_cam,vmin=my_vmin,vmax=my_vmax,interpolation='bicubic')#or bilinear, nearest, bicubic; plt.get_cmap('seismic')
plt.colorbar(im_cam,orientation='vertical')
plt.savefig("./cam_"+outputname)#how to save eps?


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