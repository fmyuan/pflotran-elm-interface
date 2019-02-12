import sys
from h5py import *
import numpy

filename = 'ic_datasets.h5'
h5file = File(filename,mode='w')

# size of domain
n = 8

# cell ids
iarray = numpy.zeros(n,'=i4')
for i in range(n):
  iarray[i] = i+1
dataset_name = 'Cell Ids'
h5dset = h5file.create_dataset(dataset_name, data=iarray)

# datasets
a = numpy.zeros(n,'=f8')

# liquid pressure
for i in range(n):
  a[i] = 101325. + 0.1*float(i+1)
dataset_name = 'Liquid_Pressure'
h5dset = h5file.create_dataset(dataset_name, data=a)

# gas saturation
for i in range(n):
  a[i] = float(i+1) / 10.
dataset_name = 'Gas_Saturation'
h5dset = h5file.create_dataset(dataset_name, data=a)

h5file.close()
print('done with everything')
  
