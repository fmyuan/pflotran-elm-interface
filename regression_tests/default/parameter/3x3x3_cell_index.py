import sys
import math
from h5py import *
import numpy as np

prefix = 'distance'

f = File(prefix+'_cell_indexed.h5',mode='w')

nx = 3
ny = 3
nz = 3

x = [0.5,1.5,2.5]

nxXny = nx*ny
n = nxXny*nz

distance = np.zeros(n,dtype=np.double)
count = 0
for k in range(nz):
    zz = x[k]
    for j in range(ny):
        yy = x[j]
        for i in range(nx):
            count += 1
            xx = x[i]
            distance[i+j*nx+k*nxXny] = math.sqrt(xx*xx+yy*yy+zz*zz)
f.create_dataset('Distance',data=distance)
iarray = np.arange(1,n+1,dtype=int)
f.create_dataset('Cell Ids',data=iarray)
f.close()

f = File(prefix+'_gridded.h5',mode='w')
h5grp = f.create_group('Distance')

nx = 4
ny = 4
nz = 4

x = [0.,1.,2.,3.]

h5grp.attrs['Dimension'] = np.string_('XYZ')
h5grp.attrs['Discretization'] = [1.,1.,1.]
h5grp.attrs['Origin'] = [0.,0.]

rarray = np.zeros((nx,ny,nz),'=f8')

for k in range(nz):
    zz = x[k]
    for j in range(ny):
        yy = x[j]
        for i in range(nx):
            xx = x[i]
            rarray[i][j][k] = math.sqrt(xx*xx+yy*yy+zz*zz)
h5grp.create_dataset('Data', data=rarray)
f.close()

print('Done')
