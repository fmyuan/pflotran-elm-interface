import sys
from h5py import *
import numpy

def create_group(h5file,group_name,cell_centered,step_interpolation,
                 origin,gradient,len_,nxy,wavelength,shift,amplitude):
    h5grp = h5file.create_group(group_name)
    if cell_centered:
        h5grp.attrs['Cell Centered'] = True
    if step_interpolation:
        h5grp.attrs['Interpolation Method'] = numpy.string_('STEP')
    h5grp.attrs['Dimension'] = numpy.string_('XY')

    if (cell_centered):
        dx = len_[0]/nxy[0]
        dy = len_[1]/nxy[1]
    else:
        dx = len_[0]/(nxy[0]-1)
        dy = len_[1]/(nxy[1]-1)

    h5grp.attrs['Origin'] = [origin[0],origin[1]]
    h5grp.attrs['Discretization'] = [dx,dy]

    rarray = numpy.zeros((nxy[0],nxy[1]),'=f8')

    twopi = 2.*numpy.pi
    y = origin[1]
    if cell_centered:
        y = 0.5*dy
    for j in range(nxy[1]):
        x = origin[0]
        if cell_centered:
            x = 0.5*dx
        for i in range(nxy[0]):
            rarray[i][j] = origin[2] + \
                           x*gradient[0] + \
                           y*gradient[1] + \
                           amplitude * \
                           numpy.sin((x+shift[0])/wavelength*twopi) * \
                           numpy.sin((y+shift[1])/wavelength*twopi)
            x += dx
        y += dy
    h5grp.create_dataset('Data',data=rarray) 

filename = 'strata_dataset.h5'
h5file = File(filename,'w')

cell_centered = False
step_interpolation = False

len_ = [20.,20.]
nxy = [41,41]

origin = [0.,0.,8.]
gradient = [-0.1,0.]
shift = [-10.,-10.]
wavelength = 40.
amplitude = 2.
create_group(h5file,'top_surface',cell_centered,step_interpolation,
             origin,gradient,len_,nxy,wavelength,shift,amplitude)

origin = [0.,0.,5.]
gradient = [-.05,0.2]
wavelength = 60.
shift = [-40.,2.]
amplitude = 3.
create_group(h5file,'middle_surface',cell_centered,step_interpolation,
             origin,gradient,len_,nxy,wavelength,shift,amplitude)

origin = [0.,0.,0.]
gradient = [0.,0.]
wavelength = 100.
shift = [10.,5.]
amplitude = 3.
create_group(h5file,'bottom_surface',cell_centered,step_interpolation,
             origin,gradient,len_,nxy,wavelength,shift,amplitude)

# for adhoc test
cell_centered = True
step_interpolation = True

origin = [0.,0.,7.5]
len_ = [20.,20.]
nxy = [41,41]

gradient = [-0.1,0.]
shift = [0.,0.]
wavelength = 20.
amplitude = 4.
create_group(h5file,'undulating_top_surface',cell_centered,step_interpolation,
             origin,gradient,len_,nxy,wavelength,shift,amplitude)

origin = [0.,0.,2.5]
gradient = [0.,0.2]
wavelength = 7.
shift = [10.,0.]
amplitude = 2.
create_group(h5file,'undulating_middle_surface',cell_centered,
             step_interpolation,
             origin,gradient,len_,nxy,wavelength,shift,amplitude)

h5file.close()

print('done')
