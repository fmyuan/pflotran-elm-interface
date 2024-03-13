import sys
from h5py import *
import numpy

filename = 'undulating_river_stage.h5'
h5file = File(filename,mode='w')

# 2d surface
h5grp = h5file.create_group('undulating_river_stage')

nx = 1
ny = 1
nz = 1
nt = 9

h5grp.attrs['Dimension'] = numpy.string_('XY')
h5grp.attrs['Cell Centered'] = True
h5grp.attrs['Discretization'] = [1e6,1e6]
h5grp.attrs['Origin'] = [0.,0.]
h5grp.attrs['Max Buffer Size'] = [2]
h5grp.attrs['Space Interpolation Method'] = numpy.string_('STEP')
h5grp.attrs['Time Interpolation Method'] = numpy.string_('LINEAR')
h5grp.attrs['Time Units'] = numpy.string_('y')

rarray = numpy.array([0.,0.25,0.5,0.75,1,1.25,1.5,1.75,2.])
h5dset = h5grp.create_dataset('Times', data=rarray)

rarray = numpy.zeros((nx,ny,nt),'=f8')

rarray[0][0][:] = numpy.array([80.,77.,79.,79.,80.,77.,79.,79.,80.])
h5dset = h5grp.create_dataset('Data', data=rarray)

h5file.close()

print('done')
