import sys
import math
from h5py import *
import numpy

filename = 'dataset.h5'
h5file = File(filename,mode='w')

nx = 5
ny = 1
nz = 1
n = nx*ny*nz

# cell indexed array where each cell id is calculated by i+j*nx+k*nx*ny
# Although indexing here is zero-based, PFLOTRAN will treat it one-based 
# internally.
cell_indexed_array = numpy.zeros((n),numpy.float64)

# write the cell indexed datasets

# time step 0
cell_indexed_array[0] = -1.e6
cell_indexed_array[1] = -1.e6
cell_indexed_array[2] = -1.e6
cell_indexed_array[3] = -1.e6
cell_indexed_array[4] = -1.e6

dataset_name = 'Pressure0'
h5dset = h5file.create_dataset(dataset_name, data=cell_indexed_array)

# time step 1
cell_indexed_array[0] = 1.098351E+05
cell_indexed_array[1] = 9.744246E+04
cell_indexed_array[2] = 8.472677E+04
cell_indexed_array[3] = 2.129996E+04
cell_indexed_array[4] = -7.909728E+05

dataset_name = 'Pressure1'
h5dset = h5file.create_dataset(dataset_name, data=cell_indexed_array)

# time step 2
cell_indexed_array[0] = 1.150173E+05
cell_indexed_array[1] = 1.050519E+05
cell_indexed_array[2] = 9.596069E+04
cell_indexed_array[3] = 7.916939E+04
cell_indexed_array[4] = -5.539533E+04

dataset_name = 'Pressure2'
h5dset = h5file.create_dataset(dataset_name, data=cell_indexed_array)

# time step 3
cell_indexed_array[0] = 1.159386E+05
cell_indexed_array[1] = 1.078158E+05
cell_indexed_array[2] = 9.969306E+04
cell_indexed_array[3] = 9.048505E+04
cell_indexed_array[4] = 5.261673E+04

dataset_name = 'Pressure3'
h5dset = h5file.create_dataset(dataset_name, data=cell_indexed_array)

# time step 4
cell_indexed_array[0] = 1.167586E+05
cell_indexed_array[1] = 1.102758E+05
cell_indexed_array[2] = 1.037930E+05
cell_indexed_array[3] = 9.743914E+04
cell_indexed_array[4] = 8.765043E+04

dataset_name = 'Pressure4'
h5dset = h5file.create_dataset(dataset_name, data=cell_indexed_array)

# time step 5
cell_indexed_array[0] = 1.180872E+05
cell_indexed_array[1] = 1.142615E+05
cell_indexed_array[2] = 1.104358E+05
cell_indexed_array[3] = 1.066101E+05
cell_indexed_array[4] = 1.033283E+05

dataset_name = 'Pressure5'
h5dset = h5file.create_dataset(dataset_name, data=cell_indexed_array)

cell_indexed_array[0] = 1.2e5
cell_indexed_array[1] = 1.2e5
cell_indexed_array[2] = 1.2e5
cell_indexed_array[3] = 1.2e5
cell_indexed_array[4] = 1.2e5

dataset_name = 'Pressure6'
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
  
