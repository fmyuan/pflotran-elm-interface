import sys
import os
try:
  pflotran_dir = os.environ['PFLOTRAN_DIR']
except KeyError:
  try:
    pflotran_dir = '../../'
  except KeyError:
    print('PFLOTRAN_DIR must point to PFLOTRAN installation directory and be defined in system environment variables.')
    sys.exit(1)
sys.path.append(pflotran_dir + '/src/python')
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
import pflotran as pft

observation_filename = 'inhibited_aerobic_respiration-obs-0.pft'

legend_fontsize = 'small'

def plot_aqueous(plt,filename,scale_string):
  plt.xlabel('Time [d]')
  plt.ylabel('Aqueous Concentration [M]')

  plt.yscale(scale_string)

  maxval = -1.e20
  minval = 1.e-10
  columns = [2,3,4]
  for icol in range(len(columns)):
    data = pft.Dataset(filename,1,columns[icol])
    ydata = data.get_array('y')
    maxval = max(maxval,np.amax(ydata))
    plt.plot(data.get_array('x'),data.get_array('y'),
             label=aq_labels[icol],c=aq_colors[icol])
  if scale_string == 'linear':
    plt.ylim(-.05*maxval,1.05*maxval)
    legend_loc = 'upper right'
  else:
    plt.ylim(0.5*minval,2.*maxval)
    legend_loc = 'lower right'

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
  # xx-small, x-small, small, medium, large, x-large, xx-large, 12, 14
  plt.legend(title='Aqueous',loc=legend_loc,fontsize=legend_fontsize)
  legend = plt.gca().get_legend()
  legend.get_frame().set_fill(False)
  legend.draw_frame(False)

def plot_immobile(plt,filename,scale_string):
  plt.twinx()
  plt.ylabel('Immobile Concentration [mol/m^3]',labelpad=20.0)
  plt.yscale(scale_string)
  maxval = -1.e20
  minval = 1.e-10
  columns = [5]
  for icol in range(len(columns)):
    data = pft.Dataset(filename,1,columns[icol])
    ydata = data.get_array('y')
    maxval = max(maxval,np.amax(ydata))
    plt.plot(data.get_array('x'),data.get_array('y'),
             label=im_labels[icol],ls='--',c=im_colors[icol])
  if scale_string == 'linear':
    plt.ylim(-0.05*maxval,1.05*maxval)
    legend_loc = 'center right'
  else:
    plt.ylim(0.5*minval,2.*maxval)
    legend_loc = 'lower left'

  plt.legend(title='Immobile',loc=legend_loc,fontsize=legend_fontsize)
  legend = plt.gca().get_legend()
  legend.get_frame().set_fill(False)
  legend.draw_frame(False)

aq_labels = []
aq_labels.append('O2(aq)')
aq_labels.append('CH2O(aq)')
aq_labels.append('CO2(aq)')

aq_colors = []
aq_colors.append('blue')
aq_colors.append('green')
aq_colors.append('red')

im_labels = []
im_labels.append('Biomass')

im_colors = []
im_colors.append('darkorange')

f = plt.figure(figsize=(10,6))

#linear scale
plt.title('Inhibited Aerobic Respiration Time History')
scale_string = 'linear'
plot_aqueous(plt,observation_filename,scale_string)
plot_immobile(plt,observation_filename,scale_string)

f.subplots_adjust(hspace=0.2,wspace=0.10,
                  bottom=.12,top=.92,
                  left=.1,right=.85)

plt.show()
