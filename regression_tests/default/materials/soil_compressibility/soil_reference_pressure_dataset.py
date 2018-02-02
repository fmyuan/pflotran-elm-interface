import sys
from h5py import *
import numpy

filename = 'soil_reference_pressure_dataset.h5'
h5file = File(filename,mode='w')

# size of domain
n = 100

# cell ids
iarray = numpy.zeros(n,'=i4')
for i in range(n):
  iarray[i] = i+1
dataset_name = 'Cell Ids'
h5dset = h5file.create_dataset(dataset_name, data=iarray)

# datasets
a = numpy.zeros(n,'=f8')

# Calcite volume fraction
for i in range(50,n):
  a[i] = 1.e5
dataset_name = 'soil_reference_pressure'
h5dset = h5file.create_dataset(dataset_name, data=a)

h5file.close()
print('done with everything')
  
