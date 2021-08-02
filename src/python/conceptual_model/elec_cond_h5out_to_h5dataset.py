import sys
from h5py import *
import numpy as np

pflotran_output_filename = 'pflotran_4x4x4.h5'
new_dataset_filename = 'ec_dataset.h5'
dataset_time_string = '10 y'

def time_string_to_time(string):
    w = string.split()
    time = float(w[0])
    if len(w) > 1:
        tu = w[1].strip()
        if tu == 's':
            # do nothing
            time *= 60.
        elif tu == 'm':
            time *= 60.
        elif tu == 'h':
            time *= 3600.
        elif tu == 'd':
            time *= 24.*3600.
        elif tu == 'w':
            time *= 7.*24.*3600.
        elif tu == 'y':
            time *= 365.*24.*3600.
        else:
            sys.exit('Unrecognized time units: {}'.format(tu))    
    return time

dataset_time = time_string_to_time(dataset_time_string)
f = File(pflotran_output_filename,mode='r')
list_of_contents = f.keys()
dataset_string = ''
for key in list_of_contents:
    if key.startswith('Time'):
        w = key.split()
        w2 = w[-2:]
        string = ' '.join(w2)
        if abs(time_string_to_time(string) - dataset_time) < 1.e-10:
            dataset_string = key + '/Electrical_Conductivity [S_m]'
            break
        
if dataset_string == '':
    sys.exit('Dataset at time {} not found.'.format(dataset_time))    
else:
    print('Copying dataset "{}" from {}'.format(dataset_string,pflotran_output_filename))
ec = f[dataset_string]
# order='F' converts the PFLOTRAN output in Fortran column 
# major format to 1D
print('Reshaping dataset')
ec_1d = np.reshape(ec,np.size(ec),order='F')
f.close()
f = File(new_dataset_filename,mode='w')
print('Writing dataset to {}'.format(new_dataset_filename))
f.create_dataset('electrical_conductivity',data=ec_1d)
iarray = np.arange(1,np.size(ec_1d)+1,dtype=np.int)
f.create_dataset('Cell Ids',data=iarray)
f.close()
print('Done')


