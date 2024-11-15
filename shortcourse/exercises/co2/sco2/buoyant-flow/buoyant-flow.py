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

name = 'buoyant-flow'
num_files = 5

files = pft.get_tec_filenames(name,range(num_files))
filenames = pft.get_full_paths(path,files)

f = plt.figure(figsize=(20,20))

plt.subplot(2,2,1)

f.suptitle(name,fontsize=16)
plt.ylabel('z (m)')
plt.xlabel('Gas Saturation ([])')

#i = 0
for ifile in range(len(filenames)):
    data = pft.Dataset(filenames[ifile],3,5)
    plt.plot(data.get_array('y'),data.get_array('x'),label=data.title)


plt.subplot(2,2,2)

f.suptitle(name,fontsize=16)
plt.ylabel('z (m)')
plt.xlabel('Temperature')

#i = 0
for ifile in range(len(filenames)):
    data = pft.Dataset(filenames[ifile],3,11)
    plt.plot(data.get_array('y'),data.get_array('x'),label=data.title)

plt.subplot(2,2,3)
f.suptitle(name,fontsize=16)
plt.ylabel('z (m)')
plt.xlabel('Aqueous CO2 Mass Fraction')

#i = 0
for ifile in range(len(filenames)):
  #i = i +1
  #if i == 4:
    data = pft.Dataset(filenames[ifile],3,6)
    plt.plot(data.get_array('y'),data.get_array('x'),label=data.title)

plt.subplot(2,2,4)

f.suptitle(name,fontsize=16)
plt.ylabel('z (m)')
#plt.ylabel('Pressure (MPa)')
plt.xlabel('Gas Pressure (MPa)')

#i = 0
for ifile in range(len(filenames)):
    data = pft.Dataset(filenames[ifile],3,4)
    plt.plot(data.get_array('y')/1e6,data.get_array('x'),label=data.title)

plt.legend(bbox_to_anchor=(1.05, 1),title='Time [d]')

plt.savefig(name+'.png')
