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
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import pflotran as pft

path = []
path.append('.')

name = 'hot_steam_injection'

files = pft.get_tec_filenames(name,range(5))
filenames = pft.get_full_paths(path,files)

f = plt.figure(figsize=(20,20))
plt.subplot(3,1,1)
f.suptitle(name,fontsize=16)
plt.xlabel('Distance')
plt.ylabel('Gas Saturation')

for ifile in range(len(filenames)):
  data = pft.Dataset(filenames[ifile],1,7)
  plt.plot(data.get_array('x'),data.get_array('y'),label=data.title)

plt.legend(loc=1,title='Time [y]')

plt.subplot(3,1,2)

f.suptitle(name,fontsize=16)
plt.xlabel('Distance')
plt.ylabel('Temperature (C)')

for ifile in range(len(filenames)):
  data = pft.Dataset(filenames[ifile],1,14)
  plt.plot(data.get_array('x'),data.get_array('y'),label=data.title)#-data2.get_array('y'),label=data.title)

plt.legend(loc=1,title='Time [y]')

plt.subplot(3,1,3)

f.suptitle(name,fontsize=16)
plt.xlabel('Distance')
#plt.ylabel('Temperature')
#plt.ylabel('Capillary Pressure')
plt.ylabel('Xal')

for ifile in range(len(filenames)):
  data = pft.Dataset(filenames[ifile],1,10) 
  plt.plot(data.get_array('x'),data.get_array('y'),label=data.title)

plt.legend(loc=1,title='Time [y]')

plt.savefig(name+'.png')
