import sys
from h5py import *
import numpy as np

new_dataset_filename = 'sec_dataset.h5'

nx = 5
ny = 4
nz = 6

nxXny = nx*ny
n = nxXny*nz

#Surface Electrical Conductivity
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

#Archie's Cementation Exponent --> m
archie_m = np.zeros(n,dtype=np.double)
count = 0
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            count += 1
            archie_m[i+j*nx+k*nxXny] = (1.9+0.001*count)
f.create_dataset('archie_cementation_exponent',data=archie_m)

#Archie's Saturation Exponent --> n
archie_n = np.zeros(n,dtype=np.double)
count = 0
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            count += 1
            archie_n[i+j*nx+k*nxXny] = (2.0+0.001*count)
f.create_dataset('archie_saturation_exponent',data=archie_n)

#Archie's Tortuosity Constant --> a
archie_a = np.zeros(n,dtype=np.double)
count = 0
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            count += 1
            archie_a[i+j*nx+k*nxXny] = (1.0+0.001*count)
f.create_dataset('archie_tortuosity_constant',data=archie_a)

iarray = np.arange(1,n+1,dtype=int)
f.create_dataset('Cell Ids',data=iarray)
f.close()

print('Done')


