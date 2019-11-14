import sys
import math
from h5py import *
import numpy

filename = 'dataset.h5'
h5file = File(filename,mode='w')

nx = 10
ny = 5
nz = 4
n = nx*ny*nz

# ijk indexed array
ijk_array = numpy.zeros((nx,ny,nz),numpy.float64)

# cell indexed array where each cell id is calculated by i+j*nx+k*nx*ny
# Although indexing here is zero-based, PFLOTRAN will treat it one-based 
# internally.
cell_indexed_array = numpy.zeros((n),numpy.float64)

# FOR ISOTROPIC PERMEABILITY ##################################################
# set values in ijk_array
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            ijk_array[i][j][k] = 1.e-12+1.e-15*(k+1)+1.e-18*(j+1)+1.e-21*(i+1)

# transfer value to cell indexed array
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            index = i + j*nx + k*nx*ny
            cell_indexed_array[index] = ijk_array[i][j][k]
            
# write the cell indexed dataset
dataset_name = 'bravo'
h5dset = h5file.create_dataset(dataset_name, data=cell_indexed_array)
            
# FOR ANISOTROPIC PERMEABILITY ################################################
ijk_array_x = ijk_array.copy()
ijk_array_y = ijk_array.copy()
ijk_array_z = ijk_array.copy()
# set values in ijk_array
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            ijk_array_x[i][j][k] = 2.e-12+1.e-15*(k+1)+1.e-18*(j+1)+1.e-21*(i+1)
            ijk_array_y[i][j][k] = 3.e-12+1.e-15*(k+1)+1.e-18*(j+1)+1.e-21*(i+1)
            ijk_array_z[i][j][k] = 4.e-12+1.e-15*(k+1)+1.e-18*(j+1)+1.e-21*(i+1)

cell_indexed_array_x = cell_indexed_array.copy()
cell_indexed_array_y = cell_indexed_array.copy()
cell_indexed_array_z = cell_indexed_array.copy()
# transfer value to cell indexed array
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            index = i + j*nx + k*nx*ny
            cell_indexed_array_x[index] = ijk_array_x[i][j][k]
            cell_indexed_array_y[index] = ijk_array_y[i][j][k]
            cell_indexed_array_z[index] = ijk_array_z[i][j][k]
            
# write the cell indexed dataset
dataset_name = 'delta'
h5dset = h5file.create_dataset(dataset_name, data=cell_indexed_array_x)
dataset_name = 'echo'
h5dset = h5file.create_dataset(dataset_name, data=cell_indexed_array_y)
dataset_name = 'foxtrot'
h5dset = h5file.create_dataset(dataset_name, data=cell_indexed_array_z)

# create a cell id array
cell_id_array = numpy.zeros((n),numpy.int32)  
for i in range(n):
  cell_id_array[i] = i+1 # add 1 for one-based indexing in PFLOTRAN
  
# write the cell id dataset
dataset_name = 'Cell Ids'
h5dset = h5file.create_dataset(dataset_name, data=cell_id_array)

h5file.close()
print('done creating cell indexed dataset')
  
