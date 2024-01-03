import sys
from h5py import *
import numpy as np

new_dataset_filename = 'sec_dataset.h5'

nx = 5
ny = 4
nz = 6

nxXny = nx*ny
n = nxXny*nz

ec = np.zeros(n,dtype=np.double)
ec[:] = 1.12
count = 0
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            count += 1
            ec[i+j*nx+k*nxXny] = (2.+0.001*count)*1.e-4
f = File(new_dataset_filename,mode='w')
f.create_dataset('surface_electrical_conductivity',data=ec)
iarray = np.arange(1,n+1,dtype=int)
f.create_dataset('Cell Ids',data=iarray)
f.close()

print('Done')


