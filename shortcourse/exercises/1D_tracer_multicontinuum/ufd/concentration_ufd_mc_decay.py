#Radionuclide transport with matrix diffusion and decay, no sorption
import sys
import os
try:
   pflotran_dir = os.environ['PFLOTRAN_DIR']
except KeyError:
   print('PFLOTRAN_DIR must point to PFLOTRAN installation directory and be defined in system environment variables.')
   sys.exit(1)

sys.path.append(pflotran_dir + '/src/python')

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
import pflotran as pft

path = []
path.append('.')

files = pft.get_tec_filenames('ufd_mc_decay',[1,3,5])
filenames = pft.get_full_paths(path,files)

f = plt.figure(figsize=(6,6))
plt.subplot(1,1,1)
f.suptitle("Am241 in Fracture",fontsize=16)
plt.xlabel('Z [m]')
plt.ylabel('Total [M]')


for ifile in range(len(filenames)):
  data = pft.Dataset(filenames[ifile],1,4)
  plt.plot(data.get_array('x'),data.get_array('y'),label='Am241 {} days'.format(data.title))

plt.legend(loc=1)
plt.savefig("conc_with_decay_fracture_am241.png")
plt.show()



plt.subplot(1,1,1)
plt.title("Np237 in Fracture",fontsize=16)
plt.xlabel('Z [m]')
plt.ylabel('Total [M]')
for ifile in range(len(filenames)):
  data = pft.Dataset(filenames[ifile],1,5)
  plt.plot(data.get_array('x'),data.get_array('y'),label='Np237 {} days'.format(data.title))

plt.legend(loc=4)
plt.savefig("conc_with_decay_fracture_np237.png")
plt.show()


files = pft.get_tec_filenames('ufd_mc_decay-sec-rank0-obs0',[1,3,5])
filenames = pft.get_full_paths(path,files)

plt.subplot(1,1,1)
plt.title("Concentration in Matrix at z = 1 [m]",fontsize=16)
plt.xlabel('X [m]')
plt.ylabel('Total [M]')


for ifile in range(len(filenames)):
  data = pft.Dataset(filenames[ifile],1,2)
  plt.plot(data.get_array('x'),np.flip(data.get_array('y'),0),label='Am241 {} days'.format(data.title))

plt.legend(loc=1)
plt.savefig("conc_with_decay_matrix.png")
plt.show()

