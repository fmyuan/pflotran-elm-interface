#Radionuclide transport without matrix diffusion
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

files = pft.get_tec_filenames('ufd_nomc',[1,3,5])
filenames = pft.get_full_paths(path,files)

f = plt.figure(figsize=(6,6))
plt.subplot(1,1,1)
f.suptitle("Am241 in Fracture",fontsize=16)
plt.xlabel('Z [m]')
plt.ylabel('Total [M]')


for ifile in range(len(filenames)):
  data = pft.Dataset(filenames[ifile],1,4)
  plt.plot(data.get_array('x'),data.get_array('y'),label='Am241 {} days'.format(data.title))


plt.legend(loc=3)
plt.savefig("conc_in_fracture_no_mc_am241.png")
plt.show()

