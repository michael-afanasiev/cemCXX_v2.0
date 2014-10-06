#! /Users/michaelafanasiev/anaconda/bin/python

import fileinput
import matplotlib
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from matplotlib import cm

plt.ion ()
plt.show ()

x = [] 
y = [] 
z = []

subsample = 1000
for line in fileinput.input ('./points.txt'):
  
  if fileinput.filelineno () % subsample == 0:
    fields = line.split ()
    x.append (float (fields[0]))
    y.append (float (fields[1]))
    z.append (float (fields[2]))
  
nLoop = 3
xPoly = []
yPoly = []
zPoly = []
for line in fileinput.input ('./facets.txt'):
  
  fields = line.split ()
  xPoly.append (float (fields[0]))
  yPoly.append (float (fields[1]))
  zPoly.append (float (fields[2]))    

facecolours = ['green', 'red', 'blue', 'yellow', 'black', 'orange', 'purple']
fig = plt.figure ()
ax  = fig.add_subplot (111, projection='3d')

ax.scatter (x, y, z, s=40, c='r')

beg = 0
end = 3
k = 0
for i in range (0, len (xPoly), nLoop):
  
  verts = [zip (xPoly[beg:end], yPoly[beg:end], zPoly[beg:end])]
  polygon = Poly3DCollection (verts, facecolor=facecolours[k], linewidths=0.1, alpha=0.4)
  ax.add_collection3d (polygon)
  plt.draw ()
  
  beg = end
  end = end + 3
  k   = (k + 1) % 7
  
  raw_input ()

