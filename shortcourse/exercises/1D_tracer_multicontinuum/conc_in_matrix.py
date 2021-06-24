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

files = pft.get_tec_filenames('tracer_1D_MC-sec-rank0-obs0',[1,2])
filenames = pft.get_full_paths(path,files)

f = plt.figure(figsize=(8,6))
plt.subplot(1,1,1)
f.suptitle("1D Tracer - Concentration in Matrix at z = 2 [m]",fontsize=16)
plt.xlabel('x [m]')
plt.ylabel('Concentration [M]')


for ifile in range(len(filenames)):
  data = pft.Dataset(filenames[ifile],1,2)
  plt.plot(data.get_array('x'),np.flip(data.get_array('y'),0),label=data.title)

plt.legend(loc=1)

plt.show()
