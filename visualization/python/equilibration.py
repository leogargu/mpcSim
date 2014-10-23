# Run this from the /python directory with:
# >> python equilibration.py 
#
# This generates an image, momentum.png, of the evolution of total mmentum at each timestep in the simulation



import matplotlib.pyplot as plt
import matplotlib.image as mpimg
#from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sys
from matplotlib.backends.backend_pdf import PdfPages


#define input and output directory
input_dir = './../../experiments/'
output_dir = input_dir
inputname='equilibration.dat'
outputname='momentum.pdf'

#get data for the slice to be plotted
data=np.genfromtxt(input_dir + inputname ,delimiter='\t')

#plot momentum
plt.plot(data[:,:-1]) #last column is the energy
plt.legend(["Jx","Jy","Jz"],loc=7)
plt.xlabel("timestep")
plt.ylabel("Total linear momentum J")
plt.title("Equilibration phase")


# Save it to disk, and add metadata
pdffig = PdfPages( output_dir + outputname )
plt.savefig(pdffig,format="pdf") 


#add metadata to figure
metadata = pdffig.infodict()
metadata['Title'] = 'Datafile plotted: '+ inputname
metadata['Author'] = 'Script used to plot this = '+ sys.argv[0]
metadata['Subject']= ''
metadata['Keywords']= ''
#metadata['Creator'] = 
#metadata['Producer']=
pdffig.close()



#add functionality to plot kinetic energy of the system  
#      *here*




