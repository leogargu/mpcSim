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
#cell_idx is contained in the upper right quadrant.
def sym_cells(ny,nz,cell_idx):
	assert isinstance(ny,int) and isinstance(nz,int),"ny,nz must be an integers"
	#include assertion here checking cell_idx is in the upper right quadrant?	
	nz_half=nz/2; #integer division
	
	#check input idx is in the upper right quadrant:
	assert cell_idx%nz >=nz_half and cell_idx <= (ny*nz_half+nz), "cell_idx should be in the upper right quadrant"
	
	alpha=2*(cell_idx%nz_half)+1;
	mu=ny*nz_half;
	if nz%2==0 and ny%2==0:
		return [cell_idx-alpha,cell_idx,2*mu-cell_idx-1,alpha+2*mu-cell_idx-1]
	else:
		B = ny*np.floor(nz*0.5);
		base= nz_half+nz*(cell_idx/nz); # integer division
		middle = ny*nz_half;
		centre = middle + nz_half;	
		#if on horizontal axis:	
		if cell_idx == nz*(cell_idx/nz)+nz_half: #exploiting integer division
			#print "on horizontal axis"
			dist = (centre - cell_idx )/nz;
			return [ cell_idx, centre + dist, centre + dist*nz ,centre - dist ]
		elif (cell_idx >=B and cell_idx <(B+nz)): #if on vertical axis
			#print "on vertical axis"
			dist = cell_idx - centre;
			return [cell_idx, centre - dist, centre - dist*nz, centre + dist*nz]
		else:
			c = 2*middle-cell_idx+nz-1;
			d = c + 2*(cell_idx-base);
			return [2*base-cell_idx,cell_idx,c,d]
		
    
#convert reference systems   
def from_YZ_to_yz(point,L_half):
	return [L_half - point[0], L_half + point[1] ]

#convert reference systems 
def from_yz_to_YZ(point,L_half):
	return [L_half - point[0], point[1] - L_half ]
	    
# Converts a local position vector into the local index of the cell that position belongs to
# Note: the convention for positions exactly at the boundary might be different from the C code.
# NOTE: the point coordinates must be in the local ref system yz (not YZ)
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
#with Y growing to the right and Z growing to the top
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
	#print "exit point found: ",B;
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
	assert volume>=0, 'Volume cannot be negative: V='+str(volume)
	assert volume <= a*a*a, 'Partial volume cannot be greater than cell volume'
	
	return volume

    	
    
# Returns a list with the arc segments specified as [A, B, h, v, cell_idx], being:
# A: first point (lower right)
# B: last point (left upper)
# h: value of the horizontal line (upper line) of the collision cell the segment is contained in
# v: value of the vertical line (left line) of the collision cell the segment is contained in
# cell_idx: local (slice-wise) index of the cell the segment is contained in
# Inputs: 
# R: radius of the circunference under consideration
# L, a, ny,nz: values of the slice considered
# The output of this function is used by the function that calculates the arc lengths and the interior volume.
def find_intersection_points(a,L,ny,nz,R):
	#initializing
	points = [];
	R_sq = R*R;
	L_half = L*0.5;    
	nz_half = nz*0.5;
	ny_half = ny*0.5;
	a_half = a*0.5;
	is_odd = False;

	#Handle cases nz=ny odd or even separately
	if nz%2!=0 or ny%2!=0 or L%2!=0 :
		is_odd = True;
		start_h = a_half;
		end_h = nz*a_half;
		end_v = a_half;
		dz = 0.0 ;
	else:
		start_h = a;
		end_h = nz*a_half;
		end_v = 0;
		dz = 0.5*a;
		
	start_v = a* math.floor( ( R - end_v ) / a ) + end_v;	#Y0
	if (R-end_v) % a == 0:
		start_v -= a	

	#Starting point: this point is the first intersection point only if nz,L,ny are even
	A=[ R , 0 ]; #starting point, rightmost cell in quadrant 1, in (Y,Z) coords [as opposed to the global (y,z) coords]
	aux_point = from_YZ_to_yz([ R, dz ],L_half);
	cell_idx = local_cellpos2idx(a,ny,nz,aux_point);
	#let's follow the circunference from A to B_last...
	v = start_v;
        h = start_h;
	while v >= 0:#end_v
   		#Find exit point (and direction) from current cell
   		B, crossing_left, crossing_up = exit_point(R_sq,h,v,A);
   		points.append([A,B,h,v,cell_idx])
   		#Update h or/and v according to exit direction. Diagonal crossings OK
   		if crossing_left:
   			if is_odd and v==end_v:
   				v = v - end_v;
   			else:
   				v = v - a;
   		if crossing_up:
   			h = h + a;
   		
   	 	#Find the index (within the slice) of the next cell
   	 	position_YZ = [v + a_half, h - a_half];
   	 	position_xy = from_YZ_to_yz(position_YZ,L_half);
   	 	cell_idx = local_cellpos2idx(a,ny,nz,position_xy);
   		#Save current B into future A
   		A=B
        #Out of while loop: add the weight for the last cell and its symmetric cells:
        #points.append([A,B_last,h,v,last_cell_idx])
        if is_odd: #exclude the first and last points (represent intersections with coordinate axis that do not coincide with grid)
        	return points[1:-1]

    	return points



def is_on_axis(cell_idx,ny,nz):
	assert ny==nz, "ny and nz should be equal";
	if ny%2==0 :
		return False;
	central_cell_idx = ny*nz/2;#integer division
	span = nz/2;
	vertical_axis = range(central_cell_idx-span,central_cell_idx+span+1,1);
	horizontal_axis = range(span,ny*nz-span,ny);
	if (cell_idx in vertical_axis) or (cell_idx in horizontal_axis):
		return True
	else:
		return False
	
    
# This function finds the volumes/length of the partial cells comprising the boundary of the cylinder   
# whose cross section is a circunference of radius R.
# quantity must be 'volume' or 'length'
def find_seg_quant(a,L,ny,nz,R,quantity):
	#Health check
	if quantity != 'volume' and quantity != 'length':
		print 'find_seg_quant: Quantity not recognised. Aborting...'
		sys.exit(1)
	assert ny==nz, "ny,nz must be equal\n";	
	#Initialise
	a_half = 0.5*a;
	q = np.empty(ny*nz); 
	if quantity == 'volume':#not technically true for cells permanently outside lumen, but can be corrected when plotting
		q.fill(a*a*a)
	else:
		q.fill(0)
	central_cell_idx = ny*nz/2;# integer division	
	#If R is too small and nz=ny is odd, there are no segments:
	if R <= a_half and ny%2!=0:
		if quantity == 'length':
			q[central_cell_idx] = 2*np.pi*R;
			return q;
		if quantity == 'volume':
			q[central_cell_idx] = 2*np.pi*R*a;
			return q;
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
   			if nz%2!=0:
   				if cell==central_cell_idx:
   					q[cell] = 4*quant;
   				elif is_on_axis(cell,ny,nz):
   					q[cell]=2*quant;
   				else:
   					q[cell]=quant
			else:
   				q[cell] = quant
   	# If ny=nz=odd number, 4 cells (E,N,W,S) need special treatment 	
	# All special cells are related by symmetry: their weight is the same.
   	if( ny%2 != 0):   	
   		#East cell index:
   		E = local_cellpos2idx(a,ny,nz,from_YZ_to_yz([R,0],0.5*L));
   		cells_idxs = sym_cells(ny,nz, E );
  		#exit point East cell
   		B = segments[0][0];
   		assert (B[1] - a_half)<1e-3, "Error handling E,W,S,O cells when nz=ny=odd\n";
   		C = [ R, 0 ];
   		v = a* math.floor( ( R - a_half ) / a ) + a_half;
   		if B[0]<=v:
   			v -= a; 
   		if quantity == 'volume':
			if B[1] != 0:
				atan_B = math.atan(B[0]/B[1]);
			else:
				atan_B = 0.5*math.pi;
			quant = a*a*(B[0]-v) - a*B[0]*B[1] + R*R*a*(0.5*math.pi-atan_B);
			#health check
			assert quant>=0, 'Volume cannot be negative: V='+str(quant)
			assert quant <= a*a*a, 'Partial volume cannot be greater than cell volume'
		elif quantity == 'length':
			quant = 2 * arc_length([0,0],R,C,B);
		#Save the same weight for those cells
   		for cell in cells_idxs:
   			q[cell] = quant;
   	if quantity=="length":
   		assert abs(q.sum()-2*np.pi*R)< 1e-5, "Error calculating lengths"
	return q
    
    
    
# Finds the weights to calculate radial averages. The weights are defined as the lengths of 
# the arc over the total lenngth of the circunference over which the radial average is performed. 
# The array of lengths can be obtained by calling the function find_seg_quant
def find_weights(lengths, R):
	weights = [];
	factor = 1.0/(2 * R * np.pi);
	for arc in lengths:
		weights.append(arc * factor)
	assert abs(sum(weights)-1.0)<1e-5 , "Error calculating weights"	
	return weights
	
	
