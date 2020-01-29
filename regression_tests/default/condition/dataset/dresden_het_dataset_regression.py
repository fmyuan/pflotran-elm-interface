import sys
from h5py import *
import numpy
import math

filename = 'parameters-16.h5'
h5in = File(filename,mode='r')

print("reading datasets")
#porosity = numpy.array(h5in['Porosity'])
permeability = numpy.array(h5in['Permeability'])

# scale by 60 to match ~ same tracer as 64m/yr.
scale = 60.
permeability[:] = permeability[:]*scale

def CalculateStats(dataset):
  mean = 0.
  maxval = -1.e20
  minval = 1.e20
  for i in range(dataset.size):
    value = dataset[i]
    mean += value
    maxval = max(value,maxval)
    minval = min(value,minval)
  mean /= dataset.size
  stdev = 0.
  for i in range(dataset.size):
    diff = dataset[i]-mean
    stdev += diff*diff
  stdev = math.sqrt(stdev/(dataset.size-1))
  return [mean, stdev, minval, maxval]

print(CalculateStats(permeability))
#print(CalculateStats(porosity))

h5in.close()

filename = 'perm_and_poros.h5'
h5file = File(filename,mode='w')

lx = 256.
ly = 256.
lz = 16.

nx = 3#256
ny = 3#256
nz = 3#16

#lx = 10.
#ly = 10.
#lz = 10.

#nx = 4
#ny = 4
#nz = 4

dx = lx / nx
dy = ly / ny
dz = lz / nz

# slightly larger dx to allow for node-centered
dx_data = 17.
nx_data = 16

#dx_data = 5.
#nx_data = 2

#cell_centered = True
cell_centered = False

dataset = numpy.zeros((nx_data*nx_data*nx_data),'=f8')
for k in range(nx_data):
  for j in range(nx_data):
    for i in range(nx_data):
      index = i + j*nx_data + k*nx_data*nx_data
      dataset[index] = i + j*10 + k*100

def GetIndex(x,dx):
  i = -999
  offset = 0.
  if cell_centered:
    offset = 0.5
  if x < offset*dx:
    i = -1
  elif x > (nx_data-offset)*dx:
    i = nx_data
  else:
    for ii in range(nx_data):
      if x >= (ii+offset)*dx and x <= (ii+1+offset)*dx:
        i = ii
        break 
  if i == -999:
    sys.exit("Error in GetIndex")
  return i

def SampleDataset(x,y,z,data):
  dx_local = dx_data
  dy_local = dx_data
  dz_local = dx_data
  i = GetIndex(x,dx_local)
  j = GetIndex(y,dy_local)
  k = GetIndex(z,dz_local)
#  print('before {} {} {}'.format(i,j,k))

  ip1 = min(i+1,nx_data-1)
  i = max(i,0)
  if cell_centered:
    x0 = (i+0.5)*dx_local
  else:
    x0 = i*dx_local
  x1 = x0
  if ip1 > i:
    x1 = x0 + dx_local

  jp1 = min(j+1,nx_data-1)
  j = max(j,0)
  if cell_centered:
    y0 = (j+0.5)*dy_local
  else:
    y0 = j*dy_local
  y1 = y0
  if jp1 > j:
    y1 = y0 + dy_local

  kp1 = min(k+1,nx_data-1)
  k = max(k,0)
  if cell_centered:
    z0 = (k+0.5)*dz_local
  else:
    z0 = k*dz_local
  z1 = z0
  if kp1 > k:
    z1 = z0 + dz_local

#  print('after {} {} {}'.format(i,j,k))

  index = i + j*nx_data + k*nx_data*nx_data
  c000 = data[index]
  index = ip1 + j*nx_data + k*nx_data*nx_data
  c100 = data[index]
  index = i + jp1*nx_data + k*nx_data*nx_data
  c010 = data[index]
  index = ip1 + jp1*nx_data + k*nx_data*nx_data
  c110 = data[index]
  index = i + j*nx_data + kp1*nx_data*nx_data
  c001 = data[index]
  index = ip1 + j*nx_data + kp1*nx_data*nx_data
  c101 = data[index]
  index = i + jp1*nx_data + kp1*nx_data*nx_data
  c011 = data[index]
  index = ip1 + jp1*nx_data + kp1*nx_data*nx_data
  c111 = data[index]

  if x1 > x0:
    xd = (x-x0)/(x1-x0)
    c00 = c000*(1-xd) + c100*xd
    c01 = c001*(1-xd) + c101*xd
    c10 = c010*(1-xd) + c110*xd
    c11 = c011*(1-xd) + c111*xd
  else:
    c00 = c000
    c01 = c001
    c10 = c010
    c11 = c011

  if y1 > y0:
    yd = (y-y0)/(y1-y0)
    c0 = c00*(1-yd) + c10*yd
    c1 = c01*(1-yd) + c11*yd
  else:
    c0 = c00
    c1 = c01

  if z1 > z0:
    zd = (z-z0)/(z1-z0)
    c = c0*(1-zd) + c1*zd
  else:
    c = c0

  return c

rarray = numpy.zeros((nx,ny,nz),'=f8')

print("writing permeability")
h5grp = h5file.create_group('Permeability')
h5grp.attrs['Dimension'] = numpy.string_('XYZ')
h5grp.attrs['Discretization'] = [dx,dy,dz]
h5grp.attrs['Origin'] = [0.,0.,0.]
h5grp.attrs['Interpolation Method'] = numpy.string_('Linear')
h5grp.attrs['Cell Centered'] = numpy.string_('True')

index = 0
z = 0.5*dz
for k in range(nz):
  print(k)
  y = 0.5*dy
  for j in range(ny):
    x = 0.5*dx
    for i in range(nx):
      rarray[i][j][k] = SampleDataset(x,y,z,permeability)
      index += 1
      x += dx
    y += dy
  z += dz
#print(rarray)
h5dset = h5grp.create_dataset('Data', data=rarray)

#print("writing porosity")
#h5grp = h5file.create_group('Porosity')
#h5grp.attrs['Dimension'] = numpy.string_('XYZ')
#h5grp.attrs['Discretization'] = [dx,dy,dz]
#h5grp.attrs['Origin'] = [0.,0.,0.]
#h5grp.attrs['Interpolation Method'] = numpy.string_('Linear')
#h5grp.attrs['Cell Centered'] = numpy.string_('True')
#
#index = 0
#for k in range(nz):
#  print(k)
#  for j in range(ny):
#    for i in range(nx):
#      rarray[i][j][k] = SampleDataset(x,y,z,porosity)
#      index += 1
#h5dset = h5grp.create_dataset('Data', data=rarray)

h5file.close()

print('done')
