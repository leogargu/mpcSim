# mpc-1.0 python script for paraview shell
try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

#import vessel_geometry
#CYLINDER
#Lx=2
#L=2
#a=1
#radius=0.5*(L-a)
#L_half=0.5*L

import sys
sys.path.append('/home/maths/markaw/Documents/mpcSim-1.0/DATA/')

from vessel_geometry import *

#import vessel_geometry

#print vessel_geometry.Lx
#print vessel_geometry.L
#print vessel_geometry.L_half
#print num_frames

Cylinder1 = Cylinder()
RenderView1 = GetRenderView()
DataRepresentation1 = Show()

Cylinder1.Center = [0,0,0] #[0.5*Lx, L_half, L_half]
Cylinder1.Height = Lx
Cylinder1.Resolution = 50
Cylinder1.Capping = 0
Cylinder1.Radius = radius

DataRepresentation1.Position = [0.5*Lx,L_half,L_half]  #[0.0, L, 0.0]
DataRepresentation1.Orientation = [0.0, 0.0, 270.0]


DataRepresentation1.Opacity = 0.55
Show()
Render()
#Show()


#SIMULATION BOX AND COLLISION CELLS
collision_grid_vtk = LegacyVTKReader( FileNames=['./collision_grid.vtk'] )

DataRepresentationgrid = Show()
DataRepresentationgrid.Representation = 'Wireframe'
RenderView1.CenterOfRotation = [L_half, L_half, L_half]
Show()
Render()

# PARTICLES
frames=['plasma0.vtk']

for i in range(1,num_frames):
	frames.append('plasma'+str(i)+'.vtk')
	


positions000 = LegacyVTKReader( FileNames= frames )
AnimationScene1 = GetAnimationScene()
AnimationScene1.EndTime = num_frames #9.0
AnimationScene1.PlayMode = 'Snap To TimeSteps'

my_representation17 = GetDisplayProperties(positions000)
my_representation17.DiffuseColor = [0.0, 0.0, 0.0]

Show()
Render()

# SHOW PARTICLES AS SPHERES 
#Glyph3 = Glyph( GlyphType="Sphere", GlyphTransform="Transform2" )
#Glyph3.Scalars = ['POINTS', '']
#Glyph3.SetScaleFactor = 0.1
#Glyph3.ScaleMode = 'scalar'
#Glyph3.GlyphType.Radius = 0.2
#my_representation73 = GetDisplayProperties(Glyph3)
#my_representation73.DiffuseColor = [1.0, 0.0, 0.0]

#Show()
#Render()

# VELOCITY FIELD
#SetActiveSource(positions000)
#Glyph5 = Glyph( GlyphType="Arrow", GlyphTransform="Transform2" )
#Glyph5.Vectors = ['POINTS', 'velocity']
#Glyph5.SetScaleFactor = 0.025

#Show()
#Render()
