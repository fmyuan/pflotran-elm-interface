import sys
import math
from h5py import *
import numpy

filename = 'dataset.h5'
h5file = File(filename,mode='w')

nx = 10
ny = 5
nz = 3
n = nx*ny*nz

# ijk indexed array
ijk_array = numpy.zeros((nx,ny,nz),numpy.float64)

# cell indexed array where each cell id is calculated by i+j*nx+k*nx*ny
# Although indexing here is zero-based, PFLOTRAN will treat it one-based 
# internally.
cell_indexed_array = numpy.zeros((n),numpy.float64)

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
dataset_name = 'Permeability'
h5dset = h5file.create_dataset(dataset_name, data=cell_indexed_array)
            
# create a cell id array
cell_id_array = numpy.zeros((n),numpy.int32)  
for i in range(n):
  cell_id_array[i] = i+1 # add 1 for one-based indexing in PFLOTRAN
  
# write the cell id dataset
dataset_name = 'Cell Ids'
h5dset = h5file.create_dataset(dataset_name, data=cell_id_array)

h5file.close()
print('done creating cell indexed dataset')
  
