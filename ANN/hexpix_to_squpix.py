#!/usr/bin/python

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import math

def my_circle_scatter(axes, x_array, y_array, radius=0.5, **kwargs):
    for x, y in zip(x_array, y_array):
        circle = plt.Circle((x,y), radius=radius, **kwargs)
        axes.add_patch(circle)
    return True

def my_polygon_q_scatter(axes, x_array, y_array, radius=0.5, **kwargs):
    ''' resolution is number of sides of polygon '''
    for x, y in zip(x_array, y_array):
		  polygon = patches.RegularPolygon((x,y), 4, radius=radius, **kwargs)
		  axes.add_patch(polygon)
    return True

def my_polygon_scatter(axes, x_array, y_array, radius=0.5, **kwargs):
    ''' resolution is number of sides of polygon '''
    for x, y in zip(x_array, y_array):
		  polygon = patches.RegularPolygon((x,y), 6, radius=radius, **kwargs)
		  axes.add_patch(polygon)
    return True	
	  
#1. produce the hexagonal grid of pixel centers

#square camera of side = cam_side
cam_side = 40

#1. start with a simple cartesian squared grid
#squ =np.array(np.arange(-20, 20, 1),np.arange(-20, 20, 1))
#squ = np.mgrid[-2:2,-2:2]
squ = np.arange(-2,3,1)
print(squ)

for x in squ:
	for y in squ:
		print(x,y)

#2. produce the hexagonal grid of pixel centers
vxq = np.zeros(25)
vyq = np.zeros(25)
vxh = np.zeros(25)
vyh = np.zeros(25)

i=0
for x in squ:
	for y in squ:
		vxq[i]=x
		vyq[i]=y
		vxh[i]=x
		vyh[i]=y-0.5*x
		i=i+1

#plot the hegonal pixels camera

fig, ax = plt.subplots(figsize=(6.,6.), facecolor="white" )
my_polygon_q_scatter(ax,vxq,vyq,radius=math.sqrt(2)/2, orientation=math.pi/4, alpha=0.2, color='r')
plt.scatter(vxq,vyq,c='b',s=4,marker='o',zorder=30)
plt.grid(True)
#plt.axis('scaled')

fig, ax = plt.subplots(figsize=(6.,6.), facecolor="white" )
my_polygon_scatter(ax,vxh,vyh,radius=1/math.sqrt(3), orientation=math.pi/2, alpha=0.2, color='r')
plt.scatter(vxh,vyh,c='b',s=4,marker='o',zorder=30)
plt.grid(True)
plt.axis('scaled')
plt.show() 
#3. rotate the hexagonally tesselleted camera to an equivalent square pixel camera

#plot the square pixels camera, with some additional zeros...

#
