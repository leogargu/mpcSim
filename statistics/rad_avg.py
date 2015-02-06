import radial_helpers
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
#import math
import sys
	
	
#http://en.wikipedia.org/wiki/Coefficient_of_determination
#http://terpconnect.umd.edu/~toh/spectrum/CurveFitting.html
def calculate_Rsq(yvals, residuals):
	y_mean=yvals.mean()
	SStot=0
	#SSreg=0
	SSres=0
	for i in range(0,len(yvals)):
		SStot += (yvals[i]-y_mean)**2
		#SSreg += (fit_points[i]-y_mean)**2
		SSres += (residuals[i])**2
		Rsq=1-(SSres/SStot)
	return Rsq


#xvals and yvals must be numpy 1D arrays
#This uses least squares:
#http://docs.scipy.org/doc/numpy/reference/generated/numpy.polyfit.html
def calculate_parabolic_fit(xvals,yvals):
	fit_coeffs, cov = np.polyfit(xvals, yvals, 2, full=False, cov=True)
	#create fitting polynomial
	fit = np.poly1d(fit_coeffs)
	fit_points=[]
	for x in xvals:
		fit_points.append(fit(x))
	residuals = yvals - fit_points
	R_sq = calculate_Rsq(yvals,residuals)
	return fit, R_sq
    
    
#auxiliary function to format matedata 
def sign_to_str(p):
	if p < 0 :
		return '-'
	else:
		return '+'
	
##############################################
# Main program
##############################################

# check number of calling arguments
if len(sys.argv)<6:
	print "call as: \n>>python rad_avg.py <dir, relative to calling directory> <name of data file> <Value of a in the simulation> <Value of L in the simulation> <Number of circunferences over which to perform the radial average><(Optional) Name of the output figure>\n"
	sys.exit(1)


#define input and output directory
input_dir = sys.argv[1]# './../experiments/'
output_dir = input_dir

#Get arguments
input_file = sys.argv[2];
a = float(sys.argv[3]);
L = float(sys.argv[4]);
R = (L-a)*0.5;
L_half = 0.5*L;
num_points = int(sys.argv[5])/2;
if len(sys.argv)==7:
	outputname=sys.argv[6]+".pdf" #change the extension if required
	outputname_density=sys.argv[6]+"_density.pdf"
else:
	outputname="output.pdf"
	outputname_density = 'output_density.pdf'
	print "Plots will be saved as output.pdf and output_density.pdf\n"

#Note that the number of points is not strictly followed (+-1)

#get header
fh = open( input_dir + input_file,'r')
header=fh.readline()
fh.close()

# process dimensions of simulation box
dimensions=header.split('\t')[:8]
nx=int(dimensions[0])
ny=int(dimensions[1])
nz=int(dimensions[2])

first_file = int(dimensions[5])		#first datafile number used in the average that generated the datafile plotted here 
stride = int(dimensions[6])		#stride used when calculating the average datafile plotted here
last_file = int(dimensions[7])	
	
num_samples = ((last_file-first_file)/stride)+1;
assert isinstance(num_samples,int), "rad_avg.py.py: error calculating number of samples"

#get slice data
data = pd.read_table(input_dir+input_file,sep='\t',skiprows=1,skipinitialspace=True,header=None);
#for diagnostic purposes mainly, we'll also do the radial average of the particle density
num_particles = data.values[:,0];

#Get radial average points
y_values = []
x_vals = []
for i in range(1,num_points+1):
	radius =  i*R/(1.0*num_points); #a + i*(R-a)/(1.0*num_points-1.0);
	x_vals.insert(0,L_half-radius)
	x_vals.append(L_half+radius)
	lengths = radial_helpers.find_seg_quant(a,L,ny,nz,radius,'length')
	warray=radial_helpers.find_weights(lengths,radius)
	y_values.append(np.dot(np.nan_to_num(data.values[:,1]),warray))
y_vals=np.concatenate((y_values[::-1],y_values))


#Calculate the best parabolic fit (least squares)
fit,R2 = calculate_parabolic_fit(x_vals,y_vals)
fit_points=[]
fit_x_values=np.linspace(0 ,ny , 100)
for x in fit_x_values:
    fit_points.append(fit(x))    
    
   
print "R^2 = ", R2
print "Slip (from data) = ", np.array(y_vals).min()
print "Slip (from fit) =", fit(0.5*a)
print "Fit: p[1]="  + str(fit[0]) + ", p[x]= "+ str(fit[1])+', p[x^2]= '+ str(fit[2])

#For latex:
print num_samples , " & " , '%.6f'%fit[0] , " & ", '%.6f'%fit[1], ' & ', '%.6f'%fit[2], " & " , '%.6f'%R2 , " \\\\"
    
#Plotting
plt.figure(1)
plt.plot(x_vals,y_vals, 'ro', markersize=3)
fig=plt.gcf()
ax=fig.gca()
ax.plot(fit_x_values,fit_points)

y1=y_vals.max()
ax.set_ylim([0,1.1*y1])
plt.vlines(range(0,ny,1),0,1.1*y1)
plt.vlines([0.5,ny-0.5],0,1.1*y1,'g')
plt.vlines([ny*0.5],0,1.1*y1,'g',linestyle='dashed');
ax.set_xlim([0,L]);
ax.set_xlim(ax.get_xlim()[::-1]) 


ax.set_xlabel('Y (cells)',fontweight="bold")
ax.set_ylabel('v_x',fontweight="bold") #Latex?		


#save averages plot to disk, prepare metadata object first
pdffig = PdfPages( output_dir + outputname )
plt.savefig( pdffig, format="pdf")#output_dir + outputname )

#add metadata to figure
metadata = pdffig.infodict()
metadata['Title'] = 'Data plotted =' + input_dir + input_file #input_file 
metadata['Author'] = 'Script used to plot this = rad_avg.py.py'


metadata['Subject']= ' Radial averages and least squares parabolic fit = '+ str(fit[0]) + sign_to_str(fit[1]) +str(abs(fit[1]))+'x'+ sign_to_str(fit[2]) + str(abs(fit[2]))+'x^2'
metadata['Keywords']= 'R^2 = '+ str(R2)+' , Slip (fit)= '+str(fit(0.5*a))+', Slip (data)= '+str(np.array(y_vals).min())
#metadata['Creator'] = 
#metadata['Producer']=


#Save figure
pdffig.close()


##########################################
#Now the particle densities
plt.figure(2)

volumes = radial_helpers.find_seg_quant(a,L,ny,nz,R,'volume')	
num_particles_cell = np.divide(num_particles,num_samples)
particle_densities = np.divide(num_particles_cell, volumes)
	
	
#Get radial average points
y_values = []
x_vals = []
num_points *= 3
a_half=0.5*a;
for i in range(1,num_points+1):
	radius =  i*R/(1.0*num_points);#a_half + i*(R-a_half)/(1.0*num_points);#i*R/(1.0*num_points);
	x_vals.insert(0,L_half-radius)
	x_vals.append(L_half+radius)
	lengths = radial_helpers.find_seg_quant(a,L,ny,nz,radius,'length')
	warray = radial_helpers.find_weights(lengths,radius)
	y_values.append(np.dot(particle_densities,warray))
y_vals=np.concatenate((y_values[::-1],y_values))


plt.plot(x_vals,y_vals, 'ro', markersize=3)
fig=plt.gcf()
ax=fig.gca()

yaux=y_vals.mean()
y1=yaux*1.8

#y1=y_vals.max()
ax.set_ylim([0,1.1*y1])
plt.vlines(range(0,ny,1),0,1.1*y1)
plt.vlines([0.5,ny-0.5],0,1.1*y1,'g')
plt.vlines([ny*0.5],0,1.1*y1,'g',linestyle='dashed');
ax.set_xlim([0,L])
ax.set_xlim(ax.get_xlim()[::-1]) 


ax.set_xlabel('Y (cells)',fontweight="bold")
ax.set_ylabel('particle density',fontweight="bold") #Latex?		



#save averages plot to disk, prepare metadata object first
pdffig = PdfPages( output_dir + outputname_density )
plt.savefig( pdffig, format="pdf")

#add metadata to figure
metadata = pdffig.infodict()
metadata['Title'] = 'Data plotted =' + input_dir + input_file #input_file 
metadata['Author'] = 'Script used to plot this = rad_avg.py.py and radial_helpers.py'


metadata['Subject']= ' Radial averages of the particle density '
#metadata['Keywords']= ''
#metadata['Creator'] = 
#metadata['Producer']=


#Save figure
pdffig.close()











