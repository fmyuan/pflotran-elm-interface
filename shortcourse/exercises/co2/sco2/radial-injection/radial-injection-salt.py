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

name = 'radial-injection-salt'
num_files = 6

files = pft.get_tec_filenames(name,range(num_files))
filenames = pft.get_full_paths(path,files)

f = plt.figure(figsize=(20,20))

plt.subplot(2,2,1)

f.suptitle(name,fontsize=16)
plt.xlabel('r (m)')
plt.ylabel('Gas Saturation ([])')

#i = 0
for ifile in range(len(filenames)):
    data = pft.Dataset(filenames[ifile],1,5)
    plt.semilogx(data.get_array('x'),data.get_array('y'),label=data.title)


plt.subplot(2,2,2)

f.suptitle(name,fontsize=16)
plt.xlabel('r (m)')
plt.ylabel('Precipitate Saturation ([])')

#i = 0
for ifile in range(len(filenames)):
    data = pft.Dataset(filenames[ifile],1,6)
    plt.semilogx(data.get_array('x'),data.get_array('y'),label=data.title)

plt.subplot(2,2,3)
f.suptitle(name,fontsize=16)
plt.xlabel('r (m)')
plt.ylabel('Aqueous CO2 Mass Fraction')

#i = 0
for ifile in range(len(filenames)):
  #i = i +1
  #if i == 4:
    data = pft.Dataset(filenames[ifile],1,7)
    plt.semilogx(data.get_array('x'),data.get_array('y'),label=data.title)

plt.subplot(2,2,4)

f.suptitle(name,fontsize=16)
plt.xlabel('r (m)')
plt.ylabel('Temperature(C)')

#i = 0
for ifile in range(len(filenames)):
    data = pft.Dataset(filenames[ifile],1,4)
    plt.semilogx(data.get_array('x'),data.get_array('y'),label=data.title)

plt.legend(bbox_to_anchor=(1.05, 1),title='Time [yr]')

plt.savefig(name+'.png')
