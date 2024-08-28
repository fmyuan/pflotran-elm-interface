import sys
from h5py import *
import numpy

def createGroup(h5file,cell_centered,step_interpolation,
                group_name,strings,dimension,discretization,origin,dataset):
  h5grp = h5file.create_group(group_name+''.join(strings))
  if cell_centered:
    h5grp.attrs['Cell Centered'] = True
  if step_interpolation:
    h5grp.attrs['Space Interpolation Method'] = numpy.string_('STEP')
  h5grp.attrs['Dimension'] = numpy.string_(dimension)
  h5grp.attrs['Discretization'] = discretization
  h5grp.attrs['Origin'] = origin
  h5dset = h5grp.create_dataset('Data', data=dataset)

def createSet(cell_centered,step_interpolation):
  nX = 11
  offset = 0.
  strings = []
  if cell_centered:
    nX -= 1
    offset = 0.5
    strings.append('_cell')
  else:
    strings.append('_node')
  if step_interpolation:
    strings.append('_step')
  else:
    strings.append('_linear')

  p_ref = 0.
  grid_spacing = 1.
  origin1 = [0.]
  origin2 = [0.,0.]
  origin3 = [0.,0.,0.]
  discretization1 = [grid_spacing]
  discretization2 = [grid_spacing,grid_spacing]
  discretization3 = [grid_spacing,grid_spacing,grid_spacing]
  rarray1 = numpy.zeros(nX,'=f8')
  rarray2 = numpy.zeros((nX,nX),'=f8')
  rarray3 = numpy.zeros((nX,nX,nX),'=f8')
  for k in range(nX):
    for j in range(nX):
      for i in range(nX):
        rarray3[i][j][k] = p_ref + (i+offset)*grid_spacing + \
                                   (j+offset)*grid_spacing*100. + \
                                   (k+offset)*grid_spacing*10000.

  createGroup(h5file,cell_centered,step_interpolation,
              'xyz',strings,'XYZ',discretization3,origin3,rarray3)
  
  k = 0
  for j in range(nX):
    for i in range(nX):
      rarray2[i][j] = rarray3[i][j][k]
  createGroup(h5file,cell_centered,step_interpolation,
              'xy',strings,'XY',discretization2,origin2,rarray2)
  
  j = 0
  for k in range(nX):
    for i in range(nX):
      rarray2[i][k] = rarray3[i][j][k]
  createGroup(h5file,cell_centered,step_interpolation,
              'xz',strings,'XZ',discretization2,origin2,rarray2)
  
  i = 0
  for k in range(nX):
    for j in range(nX):
      rarray2[j][k] = rarray3[i][j][k]
  createGroup(h5file,cell_centered,step_interpolation,
              'yz',strings,'YZ',discretization2,origin2,rarray2)
  
  j = 0
  k = 0
  for i in range(nX):
    rarray1[i] = rarray3[i][j][k]
  createGroup(h5file,cell_centered,step_interpolation,
              'x',strings,'X',discretization1,origin1,rarray1)
  
  i = 0
  k = 0
  for j in range(nX):
    rarray1[j] = rarray3[i][j][k]
  createGroup(h5file,cell_centered,step_interpolation,
              'y',strings,'Y',discretization1,origin1,rarray1)
  
  i = 0
  j = 0
  for k in range(nX):
    rarray1[k] = rarray3[i][j][k]
  createGroup(h5file,cell_centered,step_interpolation,
              'z',strings,'Z',discretization1,origin1,rarray1)

filename = 'test_datasets.h5'
h5file = File(filename,mode='w')
createSet(True,True)
createSet(False,True)
createSet(True,False)
createSet(False,False)
h5file.close()

print('done')
