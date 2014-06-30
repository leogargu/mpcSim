import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math
import sys

##############################################
#Auxiliary functions
##############################################
#Finds the local indices of the cells symmetric to the one given by the local index call_idx
def sym_cells(ny,nz,cell_idx):
	if nz%2!=0 or ny%2!=0 :
		print "nz, ny must be even. Aborting...\n";
		sys.exit(1);
	nz_half=nz/2;
	alpha=2*(cell_idx%nz_half)+1;
	mu=ny*nz_half;
	
	return [cell_idx-alpha,cell_idx,2*mu-cell_idx-1,alpha+2*mu-cell_idx-1]
    
    
# Converts a local position vector into the local index of the cell that position belongs to
# Note: the convention for positions exactly at teh boundary might be different from the C code.
def local_cellpos2idx(a,ny,nz,point):
	ainv=1.0/a;
	coords=[int(point[0]*ainv),int(point[1]*ainv)];
	
	return coords[0]*nz+coords[1]
    
    
# Calculates the arc length between two points p1 and p2 on the same circunference of radius r
# and center c (the centre is needed only to check that the input points are indeed on the circunference)
def arc_length(c,r,p1,p2):
	distSQ1=(c[0]-p1[0])*(c[0]-p1[0])+(c[1]-p1[1])*(c[1]-p1[1]);
	distSQ2=(c[0]-p2[0])*(c[0]-p2[0])+(c[1]-p2[1])*(c[1]-p2[1]);
	if (distSQ1-distSQ2) >1e-4:
		print "arc_length: p1 and p2 are not on the same circunference.Aborting...\n"
		sys.exit(1);
	norm=math.sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1]));
	rinv=1.0/r;
	
	return 2*r*math.asin(0.5*norm*rinv)



#There is A LOT of repetition in this code: clean it up.
def find_weights(a,L,ny,nz,radius):
	if nz%2!=0 or ny%2!=0 or L%2!=0 :
		print "nz, ny must be even. Aborting...\n";
		sys.exit(1);
	weights=np.empty(ny*nz); 
	weights.fill(0);    
	L_half=L/2;    
	nz_half=nz/2;
	p1=[L_half-radius,L_half]; #starting point, rightmost cell in quadrant 1
	plast=[L_half,L_half+radius];
	cell_last=local_cellpos2idx(a,ny,nz,[plast[0]-0.1*a,plast[1]-0.1*a]);
	cell_idx= local_cellpos2idx(a,ny,nz,p1);
	center=[L_half,L_half];
    
  	h=L_half+a;
  	v=int(p1[0]/a)+a;
    
   	#points=[p1];
    
   	while cell_idx != cell_last:
   		#print cell_idx
   	 	crossing_left=0;
   	 	#exiting current cell from above?
   	 	yaux=radius*radius-(h-L_half)*(h-L_half);
   	 	if yaux>=0:
   	 		y= L_half - math.sqrt( yaux ) ;
   	 		if abs(y-v)==0:
   	 			#print "diagonal crossing!\n"
   	 	 	 	p2=[v,h];
   	 	 	 	weight=arc_length(center,radius,p1,p2)/(2*radius*math.pi);
   	 	 	 	cells=sym_cells(ny,nz,cell_idx);
   	 	 	 	for cell in cells:
   	 	 	 		weights[cell]=weight;
   	 	 	 	p1=p2;
   	 	 	 	cell_idx=cell_idx+nz+1;
   	 	 	 	h=h+a;
   	 	 	 	v=v+a;
   	 	 	 	#points.append(p1);
   	 	 	 	continue;
   	 	 	if y<v and y>=(v-a):#...then I am exiting the current cell from above
   	 	 		p2=[y,h];
   	 	 	 	#calculate the arc length between p1 and p2
   	 	 	 	weight=arc_length(center,radius,p1,p2)/(2*radius*math.pi);
   	 	 	 	cells=sym_cells(ny,nz,cell_idx);
   	 	 	 	for cell in cells:
   	 	 	 		weights[cell]=weight;
   	 	 	 	cell_idx=cell_idx+1;
   	 	 	 	p1=p2; #python dangerous?
   	 	 	 	h=h+a;
   	 	 	else:
   	 	 		crossing_left=1;
   	 	else:
   	 		crossing_left=1;
        
        
		if crossing_left: #I am exiting the current cell from the left
			zaux= radius*radius - (v-L_half)*(v-L_half)
			if zaux>=0:
				z= L_half + math.sqrt( zaux );
				if abs(z-h)==0:
					#print "diagonal crossing!\n"
					p2=[v,h];
					weight=arc_length(center,radius,p1,p2)/(2*radius*math.pi);
					cells=sym_cells(ny,nz,cell_idx);
					for cell in cells:
						weights[cell]=weight;
					p1=p2;
					cell_idx=cell_idx+nz+1;
					h=h+a;
					v=v+a;
					#points.append(p1);
					continue;
                    
        	        	if z<=h and z>(h-a):
        	        		p2=[v,z];
        	        		#calculate arc length between p1 and p2
        	        		weight=arc_length(center,radius,p1,p2)/(2*radius*math.pi);
        	        		cells=sym_cells(ny,nz,cell_idx);
        	        		for cell in cells:
        	        			weights[cell]=weight;
        	        		cell_idx=cell_idx+nz;
        	        		p1=p2;
        	        		v=v+a;
        	        	else:
        	        		print "Error. Attempted to cross left but couldn't. Aborting...\n"
        	        		sys.exit(1);
        	        else:
        	        	print "Error. Not knowing where to cross. Aborting...\n"
        	        	sys.exit(1);
        	#points.append(p1);
        #Out of while loop: add the arc in the last cell:
        weight=arc_length(center,radius,p1,plast)/(2*radius*math.pi);
        cells=sym_cells(ny,nz,cell_idx);
        for cell in cells:
        	weights[cell]=weight;
        #points.append(plast);
        
    	#print cell_idx
    
    	return weights #,points
    
##############################################
# Main program
##############################################

# check number of calling arguments
if len(sys.argv)<6:
	print "call as: \n>>python rad_avg.py <full path and name of data file in /experiments> <Value of a in the simulation> <Value of L in the simulation> <radius over which to perform the average><viscosity /mu/ >\n"
	sys.exit(1)


#define input and output directory
input_dir = './../experiments/'
output_dir = input_dir

#Get arguments
input_file = sys.argv[1];
a = float(sys.argv[2]);
L = float(sys.argv[3]);
radius = float(sys.argv[4]);
viscosity=float(sys.arg[5]);

#get header
fh = open( input_dir + input_file,'r')
header=fh.readline()
fh.close()

# process dimensions of simulation box
dimensions=header.split('\t')[:3]
nx=int(dimensions[0])
ny=int(dimensions[1])
nz=int(dimensions[2])


#get slice data
data = pd.read_table(input_dir+input_file,sep='\t',skiprows=1,skipinitialspace=True,header=None);

#Get weights
warray=find_weights(a,L,ny,nz,radius);

#Check weights have been calculated correctly
assert (warray.sum()-1.0) <1e-6

#Calculate average over the circunference, and print to screen
#print np.dot(np.nan_to_num(data.values[:,1]),warray)
average = np.average(np.nan_to_num(data.values[:,1]),weights=warray)
print "Weighted mean:", average;

#Calculate the weighted standard deviation
#Not clear what expression to follow, undecided about how to correct bias for small samples. For example:
#http://www.itl.nist.gov/div898/software/dataplot/refman2/ch2/weightsd.pdf
# Here will follow GNU GSL:
#http://www.gnu.org/software/gsl/manual/html_node/Weighted-Samples.html
V2=np.dot(warray,warray);
V1=warray.sum();
average_array=np.empty(ny*nz); 
average_array.fill(average); 
diffSQ=np.square( np.subtract(np.nan_to_num(data.values[:,1]) , average_array) );
print "Weighted std: ", ( V1 / (V1*V1-V2) )* np.dot(warray,diffSQ)




################################################################
#Merge me properly later - This is work in progress

#generate average points:
app_points = np.empty(10); 
app_points.fill(0);

R = (L-a)*0.5

for i in range(0,11):
	r = a + i*(R-3*a)*0.05
	warray=find_weights(a,L,ny,nz,r)
	assert (warray.sum()-1.0) <1e-6
	average = np.average(np.nan_to_num(data.values[:,1]),weights=warray)
	app_points[i]=average

Vmax=np.nanmax(data.values[:,1])
# plot several radial averages and compare with pot of videal=Vmax(1-(r/R)^2)

x=np.linspace(-R-1,R+1)
y=Vmax*(1-(x/R)*(x/R))
plt.plot(x, )


