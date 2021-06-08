# For 2D Vadose exercise: 
# Make an initial condition dataset for a model with regional sloping water table
# from previous runs to steady-state that specify the east and west boundaries.

import h5py as h5
import numpy as np


f = h5.File('./ic_dataset.h5','w')
west_ic = h5.File('./west.h5', 'r')
east_ic = h5.File('./east.h5', 'r')

#Extract variables that will be used as boundary and initial conditions

time_slice_string = 'Time:  2.20000E+02 y'
t_w = west_ic[time_slice_string]['Temperature [C]'][0,0,:]
lp_w = west_ic[time_slice_string]['Liquid_Pressure [Pa]'][0,0,:]
gp_w = west_ic[time_slice_string]['Gas_Pressure [Pa]'][0,0,:]
gs_w = west_ic[time_slice_string]['Gas_Saturation'][0,0,:]
t_e = east_ic[time_slice_string]['Temperature [C]'][0,0,:]
lp_e = east_ic[time_slice_string]['Liquid_Pressure [Pa]'][0,0,:]
gp_e = east_ic[time_slice_string]['Gas_Pressure [Pa]'][0,0,:]
gs_e = east_ic[time_slice_string]['Gas_Saturation'][0,0,:]
xgl_w = west_ic[time_slice_string]['X_g^l'][0,0,:]
xgl_e = east_ic[time_slice_string]['X_g^l'][0,0,:]

f.create_group('bounding_gas_pressure')
f.create_group('bounding_liquid_pressure')
f.create_group('bounding_temperature')
f.create_group('bounding_gas_saturation')
f.create_group('bounding_mole_fraction')

#group west and east faces in the dataset
lp = np.vstack((lp_w,lp_e))
gp = np.vstack((gp_w,gp_e))
gs = np.vstack((gs_w,gs_e))
t = np.vstack((t_w,t_e))
xgl = np.vstack((xgl_w,xgl_e))

#require non-zero gas pressure
gp = np.maximum(gp,lp)

f['bounding_liquid_pressure'].create_dataset('Data',data=lp)
f['bounding_gas_saturation'].create_dataset('Data',data=gs)
f['bounding_gas_pressure'].create_dataset('Data',data=gp)
f['bounding_temperature'].create_dataset('Data',data=t)
f['bounding_mole_fraction'].create_dataset('Data',data=xgl)

# Assign group attributes (cell-centered)
f['bounding_temperature'].attrs.create('Dimension',['XZ'],dtype='|S2')
f['bounding_temperature'].attrs.create('Discretization',[100.0,0.1])
#f['bounding_temperature'].attrs.create('Origin',[0.,0.0])
f['bounding_temperature'].attrs.create('Origin',[-50.,0.0])
f['bounding_temperature'].attrs.create('Cell Centered',[True])

f['bounding_liquid_pressure'].attrs.create('Dimension',['XZ'],dtype='|S2')
f['bounding_liquid_pressure'].attrs.create('Discretization',[100.0,0.1])
#f['bounding_temperature'].attrs.create('Origin',[0.,0.0])
f['bounding_liquid_pressure'].attrs.create('Origin',[-50.,0.0])
f['bounding_liquid_pressure'].attrs.create('Cell Centered',[True])

f['bounding_gas_pressure'].attrs.create('Dimension',['XZ'],dtype='|S2')
f['bounding_gas_pressure'].attrs.create('Discretization',[100.0,0.1])
#f['bounding_temperature'].attrs.create('Origin',[0.,0.0])
f['bounding_gas_pressure'].attrs.create('Origin',[-50.,0.0])
f['bounding_gas_pressure'].attrs.create('Cell Centered',[True])

f['bounding_gas_saturation'].attrs.create('Dimension',['XZ'],dtype='|S2')
f['bounding_gas_saturation'].attrs.create('Discretization',[100.0,0.1])
#f['bounding_temperature'].attrs.create('Origin',[0.,0.0])
f['bounding_gas_saturation'].attrs.create('Origin',[-50.,0.0])
f['bounding_gas_saturation'].attrs.create('Cell Centered', [True])

f['bounding_mole_fraction'].attrs.create('Dimension',['XZ'],dtype='|S2')
f['bounding_mole_fraction'].attrs.create('Discretization',[100.0,0.1])
#f['bounding_temperature'].attrs.create('Origin',[0.,0.0])
f['bounding_mole_fraction'].attrs.create('Origin',[-50.,0.0])
f['bounding_mole_fraction'].attrs.create('Cell Centered', [True])

f.close()
west_ic.close()
east_ic.close()

