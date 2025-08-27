import h5py
import numpy as np

'''
Purpose:    Create HDF5-format geomechanics datasets, for use in spatially-varying
            geomechanical properties.
Author:     Kyle Mosley (kyle.mosley@wsp.com)
Date:       08-2025
'''

# delete datasets in existing file
h5path = 'geomech_prop_1.h5'

# generate some (randomly-distributed) geomech properties
rd_data = 1100 + np.array(range(1, 28))
ym_data = 1000000100 + np.array(range(1, 28))
pr_data = 0.45 + np.array(range(1, 28)) * 0.0001
bc_data = 0 + np.array(range(1, 28)) * 0.0001
tec_data = 0 + np.array(range(1, 28)) * 0.0001

# add new datasets to file
with h5py.File(h5path, 'w') as fh:
    fh.create_dataset('Cell Ids', (27,), data=list(range(1, 28)))
    fh.create_dataset('Rock_density_1', (27,), data=rd_data)
    fh.create_dataset('Youngs_modulus_1', (27,), data=ym_data)
    fh.create_dataset('Poissons_ratio_1', (27,), data=pr_data)
    fh.create_dataset('Biot_coeff_1', (27,), data=bc_data)
    fh.create_dataset('Thermal_exp_coeff_1', (27,), data=tec_data)
    
    
    

# delete datasets in existing file
h5path = 'geomech_prop_2.h5'

# generate some (randomly-distributed) geomech properties
rd_data = 1200 + np.array(range(1, 28))
ym_data = 1000000200 + np.array(range(1, 28))
pr_data = 0.45 + np.array(range(1, 28)) * 0.00001
bc_data = 0 + np.array(range(1, 28)) * 0.00001
tec_data = 0 + np.array(range(1, 28)) * 0.00001

# add new datasets to file
with h5py.File(h5path, 'w') as fh:
    fh.create_dataset('Cell Ids', (27,), data=list(range(1, 28)))
    fh.create_dataset('Rock_density_2', (27,), data=rd_data)
    fh.create_dataset('Youngs_modulus_2', (27,), data=ym_data)
    fh.create_dataset('Poissons_ratio_2', (27,), data=pr_data)
    fh.create_dataset('Biot_coeff_2', (27,), data=bc_data)
    fh.create_dataset('Thermal_exp_coeff_2', (27,), data=tec_data)