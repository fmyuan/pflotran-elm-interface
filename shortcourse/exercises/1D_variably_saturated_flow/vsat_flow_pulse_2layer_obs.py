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

files = []
files.append('vsat_flow_pulse_2layer-obs-0.tec')
filenames = pft.get_full_paths(path,files)

f = plt.figure(figsize=(6,6))
plt.subplot(1,1,1)
f.suptitle("1D Variably Saturated Flow - Pulse, 2 Layer, Observation",fontsize=14)
plt.xlabel('Time [y]')
plt.ylabel('Saturation [-]')

#plt.xlim(0.,1.)
#plt.ylim(0.,1.)
#plt.grid(True)

data = pft.Dataset(filenames[0],1,3)
plt.plot(data.get_array('x'),data.get_array('y'),label='(0.5 0.5 4.95)')

#'best'         : 0, (only implemented for axis legends)
#'upper right'  : 1,
#'upper left'   : 2,
#'lower left'   : 3,
#'lower right'  : 4,
#'right'        : 5,
#'center left'  : 6,
#'center right' : 7,
#'lower center' : 8,
#'upper center' : 9,
#'center'       : 10,
plt.legend(loc=1,title='Location')
# xx-small, x-small, small, medium, large, x-large, xx-large, 12, 14
plt.setp(plt.gca().get_legend().get_texts(),fontsize='small')
#plt.setp(plt.gca().get_legend().get_texts(),linespacing=0.)
plt.gca().get_legend().get_frame().set_fill(False)
plt.gca().get_legend().draw_frame(False)
#plt.gca().yaxis.get_major_formatter().set_powerlimits((-1,1))

f.subplots_adjust(hspace=0.2,wspace=0.2,
                  bottom=.12,top=.9,
                  left=.14,right=.9)

plt.show()
