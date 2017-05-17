import sys
from h5py import *
import numpy as np

filename = 'dataset.h5'
h5file = File(filename,mode='w')

L = 1.0     # [m]
p0 = 1.0e6  # [Pa]
p_offset = 1.0e6  # [Pa]

# 1d line in x
# Pressure boundary condition; SOUTH
h5grp = h5file.create_group('x_line_node_centered_south')
h5grp.attrs['Dimension'] = np.string_('X')
# Delta length between points [m]
h5grp.attrs['Discretization'] = [0.1]
# Location of origin
h5grp.attrs['Origin'] = [0.]
# Load the dataset values
nx = 21  # nx+1
rarray = np.zeros(nx,'=f8')
for i in range(nx):
  x = (float(i)*2*L)/(nx-1)
  rarray[i] = p0*x/L + p_offset
h5dset = h5grp.create_dataset('Data', data=rarray)

# 1d line in x
# Pressure boundary condition; NORTH
h5grp = h5file.create_group('x_line_node_centered_north')
h5grp.attrs['Dimension'] = np.string_('X')
# Delta length between points [m]
h5grp.attrs['Discretization'] = [0.1]
# Location of origin
h5grp.attrs['Origin'] = [0.]
# Load the dataset values
nx = 21  # nx+1
rarray = np.zeros(nx,'=f8')
for i in range(nx):
  x = (float(i)*2*L)/(nx-1)
  rarray[i] = (p0/L)*(x+(2*L)) + p_offset
h5dset = h5grp.create_dataset('Data', data=rarray)

# 1d line in y
# Pressure boundary condition; EAST
h5grp = h5file.create_group('y_line_node_centered_east')
h5grp.attrs['Dimension'] = np.string_('Y')
# Delta length between points [m]
h5grp.attrs['Discretization'] = [0.1]
# Location of origin
h5grp.attrs['Origin'] = [0.]
# Load the dataset values
ny = 11  #ny+1
rarray = np.zeros(ny,'=f8')
for j in range(ny):
  y = (float(j)*L)/(ny-1)
  rarray[j] = (p0/L)*((2*L)+(2*y)) + p_offset
h5dset = h5grp.create_dataset('Data', data=rarray)
