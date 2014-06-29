import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math
import sys

##############################################
#Auxiliary functions

def sym_cells(ny,nz,cell_idx):
    if nz%2!=0 or ny%2!=0 :
        print "nz, ny must be even. Aborting...\n";
        sys.exit(1);
    nz_half=nz/2;
    alpha=2*(cell_idx%nz_half)+1;
    mu=ny*nz_half;
    
    return [cell_idx-alpha,cell_idx,2*mu-cell_idx-1,alpha+2*mu-cell_idx-1]
    
    
    
def local_cellpos2idx(ny,nz,point):
    ainv=1.0/a;
    coords=[int(point[0]*ainv),int(point[1]*ainv)];
    return coords[0]*nz+coords[1]
    
    
def arc_length(r,p1,p2):
    #check that these points are on the circunference?
    norm=math.sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1]));
    rinv=1.0/r;
    return 2*r*math.asin(0.5*norm*rinv)





#There is A LOT of repetition in this code: clean it up.
a=1.0;
def find_weights(ny,nz,L,radius):
    if nz%2!=0 or ny%2!=0 or L%2!=0 :
        print "nz, ny must be even. Aborting...\n";
        sys.exit(1);
        
    weights=np.empty(ny*nz); 
    weights.fill(0);    
    L_half=L/2;    
    nz_half=nz/2;
    p1=[L_half-radius,L_half]; #starting point, rightmost cell in quadrant 1
    plast=[L_half,L_half+radius];
    cell_last=local_cellpos2idx(ny,nz,[plast[0]-0.1*a,plast[1]-0.1*a]);
    cell_idx= local_cellpos2idx(ny,nz,p1);
    
    h=L_half+a;
    v=int(p1[0]/a)+a;
    
    points=[p1];
    
    while cell_idx != cell_last:
        print cell_idx
        crossing_left=0;
        #exiting current cell from above?
        yaux=radius*radius-(h-L_half)*(h-L_half);
        if yaux>=0:
            y= L_half - math.sqrt( yaux ) ;
            if abs(y-v)==0:
                #print "diagonal crossing!\n"
                p2=[v,h];
                weight=arc_length(radius,p1,p2)/(2*radius*math.pi);
                cells=sym_cells(ny,nz,cell_idx);
                for cell in cells:
                    weights[cell]=weight;
                p1=p2;
                cell_idx=cell_idx+nz+1;
                h=h+a;
                v=v+a;
                points.append(p1);
                continue;
                #sys.exit(1);
        
            if y<v and y>=(v-a):#then I am exiting the current cell from above
                p2=[y,h];
                #calculate the arc length between p1 and p2
                weight=arc_length(radius,p1,p2)/(2*radius*math.pi);
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
                    #print "diagonal crossing! (modify algorithm). Aborting...\n"
                    p2=[v,h];
                    weight=arc_length(radius,p1,p2)/(2*radius*math.pi);
                    cells=sym_cells(ny,nz,cell_idx);
                    for cell in cells:
                        weights[cell]=weight;
                    p1=p2;
                    cell_idx=cell_idx+nz+1;
                    h=h+a;
                    v=v+a;
                    points.append(p1);
                    continue;
                    #sys.exit(1);
                    
                if z<=h and z>(h-a):
                    p2=[v,z];
                    #calculate arc length between p1 and p2
                    weight=arc_length(radius,p1,p2)/(2*radius*math.pi);
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
        points.append(p1);
        
    #add the arc in the last cell:
    weight=arc_length(radius,p1,plast)/(2*radius*math.pi);
    cells=sym_cells(ny,nz,cell_idx);
    for cell in cells:
        weights[cell]=weight;
    points.append(plast);
        
    print cell_idx
    return weights,points
    
#############################################

# Give L and a as arguments

# check number of calling arguments
if len(sys.argv)<3:
	print "call as: \n>>python slice_compare.py <full path and name of SAM file in /experiments> <full path and name of CAM file in /experiments> <(optional) name output file, without extension>\n"
	sys.exit(1)


#define input and output directory
input_dir = './../../experiments/'
output_dir = input_dir

#get header
fh = open('./CAMtest.dat','r')
header=fh.readline()
fh.close()

# process dimensions of simulation box
dimensions=header.split('\t')[:3]
nx=int(dimensions[0])
ny=int(dimensions[1])
nz=int(dimensions[2])

input_file='./CAMtest.dat';

#get data for the slice to be plotted
data = pd.read_table(input_file,sep='\t',skiprows=1,skipinitialspace=True,header=None);



warray,mypoints=find_weights(ny,nz,ny,5*a)
warray.sum() 

np.dot(np.nan_to_num(data.values[:,1]),warray)



