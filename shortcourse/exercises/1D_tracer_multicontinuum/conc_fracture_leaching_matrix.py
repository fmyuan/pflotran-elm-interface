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
files=[]

files.append('tracer_1D_MC_leaching_test-obs-0.pft')
filenames = pft.get_full_paths(path,files)

f = plt.figure(figsize=(8,6))
plt.subplot(1,1,1)
f.suptitle("1D Tracer - Concentration in Fracture at z = 2 [m] with leaching matrix",fontsize=16)
plt.xlabel('Time [days]')
plt.ylabel('Concentration [M]')


data = pft.Dataset(filenames[0],1,2)
plt.plot(data.get_array('x'),data.get_array('y'),label=data.title)

plt.legend(loc=1)

plt.show()
