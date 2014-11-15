import pandas as pd
import matplotlib.pyplot as plt
#import matplotlib.image as mpimg
#from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math
import sys

#Everything in this file refers to the YZ reference system, centred at the centre of the circunference and with the Y axis inverted with respect 
# the global xy reference system

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
    
#convert reference systems   
def from_YZ_to_yz(point,L_half):
	return [L_half - point[0], L_half + point[1] ]

#convert reference systems 
def from_yz_to_YZ(point,L_half):
	return [L_half - point[0], point[1] - L_half ]
	    
# Converts a local position vector into the local index of the cell that position belongs to
# Note: the convention for positions exactly at the boundary might be different from the C code.
# NOTE: the point coordinates must be in the loocal ref system yz (not YZ)
def local_cellpos2idx(a,ny,nz,point):
	ainv=1.0/a;
	coords=[int(point[0]*ainv),int(point[1]*ainv)];
	
	return coords[0]*nz+coords[1]

# Calculates the arc length between two points p1 and p2 on the same circunference of radius r
# and center c (the centre is needed only to check that the input points are indeed on the circunference)
# This function works independently of the reference system chosen.
def arc_length(c,r,p1,p2):
	#health check
	distSQ1=(c[0]-p1[0])*(c[0]-p1[0])+(c[1]-p1[1])*(c[1]-p1[1]);
	distSQ2=(c[0]-p2[0])*(c[0]-p2[0])+(c[1]-p2[1])*(c[1]-p2[1]);
	if (distSQ1-distSQ2) > 1e-5:
		print "arc_length: p1 and p2 are not on the same circunference.Aborting...\n"
		sys.exit(1);
	#the actual calculation
	AB_length=math.sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1]));
	rinv=1.0/r;
	
	return 2*r*math.asin(0.5*AB_length*rinv)



#Returns boolean values for exiting directions from current cell, following the circunference (upper right quadrant only)
#R is the radius of the circunference, R_sq=R*R
#h and v define the current cell: h is the horizontal line above it and v is teh vertical line to its left
#This function uses YZ coordinates (as opposed to xy global coords), i.e. centre of reference system is at centre of the circunference, 
#with Y growing to the right and Z growing to the right
#It return the exit point YZ coordinates
#Interpreting the additional output of this function: 1,0 - exiting left; 0,1 -exiting above; 1,1-exiting through upper left diagonal; 0,0 - error
def exit_point(R_sq,h,v, A):
	crossing_left = 0;
	crossing_up = 0;
	h_sq=h*h;
	v_sq=v*v;
	
	#exiting current cell from above?
   	#we exit from above if there is an intersection between the h line and the circunference
	Y2aux=R_sq-h_sq;
   	if Y2aux>=0:
   		Yaux=math.sqrt(Y2aux);
   		#print "Yaux= ", Yaux, "h= ", h, "v= ", v
   	 	if Yaux >= v and Yaux < A[0]:
   	 		crossing_up=1;
   	 		B=[Yaux , h];
        #print "Y^2aux= ", Y2aux
   	#exiting current cell to the left?
   	#we exit the current cell to the left is there is a valid intersection between the circunference and the vertical line v
   	Z2aux = R_sq-v_sq;
   	if Z2aux>=0:
   		Zaux=math.sqrt(Z2aux);
   		#print "Zaux= ",Zaux, "h= ", h ,"v= ", v
   		if Zaux <= h and Zaux > A[1]:
   			crossing_left=1;
   			B=[ v, Zaux ];
        #print "Z^2aux= ", Z2aux
	#sanity check
	if crossing_left and crossing_up:
		assert B==[v,h]
	if not crossing_left and not crossing_up:
		print "Error in function exit_point. Aborting..."
		sys.exit(1);
	
	return B, crossing_left, crossing_up


#Calculates the partial volume inside the lumen.
# v, h define the current cell
def partial_volume(A, B, h, v, R, a):
	volume = a*( B[0] - v ) - (h-a)*( A[0]-B[0] ) + 0.5*( A[0]*A[1]-B[0]*B[1] )
	
	if A[1] != 0:
		atan_A = math.atan(A[0]/A[1]);
	else:
		atan_A = 0.5*math.pi;
	
	if B[1] != 0:
		atan_B = math.atan(B[0]/B[1]);
	else:
		atan_B = 0.5*math.pi;
		
	volume = a*( volume + 0.5*R*R*( atan_A - atan_B ) );
	
	#health check
	assert volume>=0, 'Volume cannot be negative'
	assert volume <= a*a*a, 'Partial volume cannot be greater than cell volume'
	
	return volume

    	
    
# Returns a list with the arc segments specified as [A, B, h, v, cell_idx], being:
# A: first point (lower right)
# B: last point (left upper)
# h: value of the horizontal line (upper line) of the collision cell the segment is contained in
# v: value of the vertical line (left line) of the collision cell the segment is contained in
# cell_idx: local (slie-wise) index of the cell the segment is contained in
# Inputs: 
# R: radius of the circunference under consideration
# L, a, ny,nz: values of the slice considered
# The output of this function is used by the function that calculates the arc lengths and the interior volume.
def find_intersection_points(a,L,ny,nz,R):
	#some health checks...
	if nz%2!=0 or ny%2!=0 or L%2!=0 :
		print "nz, ny must be even. Aborting...\n";
		sys.exit(1);
	#initializing
	points = [];
	R_sq = R*R;
	L_half = L*0.5;    
	nz_half = nz*0.5;
	ny_half = ny*0.5;
	a_half = a*0.5;
	#Starting point (first point A)
	A=[ R , 0 ]; #starting point, rightmost cell in quadrant 1, in (Y,Z) coords [as opposed to the global (y,z) coords]
	aux_point = from_YZ_to_yz([ R, 0.5*a ],L_half);
	cell_idx = local_cellpos2idx(a,ny,nz,aux_point);
	
	#end point (while loop stopping reference)
	B_last=[ 0, R ];
	aux_point = from_YZ_to_yz([ 0.5*a, R ],L_half);
	last_cell_idx = local_cellpos2idx(a,ny,nz,aux_point); 
	if R % a == 0:
		last_cell_idx -= 1
	
	#The current cell is below the horizontal line h, and to the right of the vertical line v
  	h = a;		 		#Z0+a
  	v= a* math.floor( R / a );	#Y0
  	if R % a == 0:
  		v -= 1
  	#v = (ny_half - 1)*a; 	#Y0
        
	#let's follow the circunference from A to B_last...
	while cell_idx != last_cell_idx:
   		#Find exit point (and direction) from current cell
   		B, crossing_left, crossing_up = exit_point(R_sq,h,v,A);
   		points.append([A,B,h,v,cell_idx])
   		#Update h or/and v according to exit direction. Diagonal crossings OK
   		if crossing_left:
   			v = v - a;
   		if crossing_up:
   			h = h + a;
   		# Debugging:
   		if crossing_up and crossing_left:
   			print "Diagonal crossing!"
   	 	#Find the index (within the slice) of the next cell
   	 	position_YZ = [v + a_half, h - a_half];
   	 	position_xy = from_YZ_to_yz(position_YZ,L_half);
   	 	cell_idx = local_cellpos2idx(a,ny,nz,position_xy);
   		#Save current B into future A
   		A=B
        #Out of while loop: add the weight for the last cell and its symmetric cells:
        points.append([A,B_last,h,v,last_cell_idx])
	
    	return points

    
# This function finds the volumes of the partial cells comprising the boundary of the cylinder   
# whose cross section is a circunference of radius R.
# quantity must be 'volume' or 'length'
def find_seg_quant(a,L,ny,nz,R,quantity):
	#Health check
	if quantity != 'volume' and quantity != 'length':
		print 'find_seg_quant: Quantity not recognised. Aborting...'
		sys.exit(1)
	#Initialise
	q = np.empty(ny*nz); 
	if quantity == 'volume':
		q.fill(a*a*a)
	else:
		q.fill(0)
	#Get segments
	segments = find_intersection_points(a,L,ny,nz,R);
	#Find volumes determined by the segments, and its symmetric counterparts	
	for segment in segments:
		if quantity == 'volume':
			quant = partial_volume(segment[0], segment[1], segment[2], segment[3], R, a);
		elif quantity == 'length':
			quant = arc_length([0,0],R,segment[0],segment[1])
		cells_idxs = sym_cells(ny,nz,segment[4]);
   		#Save the same weight for those cells
   		for cell in cells_idxs:
   			q[cell] = quant	
	return q
    
    
    
# Finds the weights to calculate radial averages. The weights are defined as the lengths of 
# the arc over the total lenngth of the circunference over which the radial average is performed. 
# The array of lengths can be obtained by calling the function find_seg_quant
def find_weights(lengths, R):
	weights = [];
	factor = 1.0/(2 * R * np.pi);
	for arc in lengths:
		weights.append(arc * factor)
	return weights
	
	
