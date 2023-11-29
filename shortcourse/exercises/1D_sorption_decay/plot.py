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

observation_filename = 'transport_sorption_decay-obs-0.pft'
plot_filename = 'transport_sorption_decay-002.tec'

legend_fontsize = 'small'

def plot_aq_conc_breakthrough(plt,filename,scale_string):
  plt.xlabel('Time [y]')
  plt.ylabel('Aqueous Concentration [M]')
  plt.yscale(scale_string)

  maxval = -1.e20
  minval = 1.e-10
  columns = [2,3,4,5,6,7]
  for icol in range(len(columns)):
    data = pft.Dataset(filename,1,columns[icol])
    ydata = data.get_array('y')
    maxval = max(maxval,np.amax(ydata))
    plt.plot(data.get_array('x'),data.get_array('y'),
             label=aq_labels[icol],c=aq_colors[icol],ls=aq_linestyles[icol])
  if scale_string == 'linear':
    plt.ylim(-0.02*maxval,1.05*maxval)
  else:
    plt.ylim(0.5*minval,2.*maxval)

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
  plt.legend(title='Aqueous',loc='center left',fontsize=legend_fontsize)
  legend = plt.gca().get_legend()
  legend.get_frame().set_fill(False)
  legend.draw_frame(False)

def plot_aq_conc_profile(plt,filename,scale_string):
  plt.xlabel('X [m]')
  plt.ylabel('Aqueous Concentration [M]')
  plt.yscale(scale_string)
  maxval = -1.e20
  minval = 1.e-10
  columns = [4,5,6,7,8,9]
  for icol in range(len(columns)):
    data = pft.Dataset(filename,1,columns[icol])
    ydata = data.get_array('y')
    maxval = max(maxval,np.amax(ydata))
    plt.plot(data.get_array('x'),data.get_array('y'),
             label=aq_labels[icol],c=aq_colors[icol],ls=aq_linestyles[icol])
  if scale_string == 'linear':
    ymin = -.02*maxval
    ymax = 1.05*maxval
  else:
    ymin = 0.5*minval
    ymax = 2.*maxval
  plt.ylim(ymin,ymax)
  plt.plot([50.,50.],[ymin,ymax],c='grey',ls='dotted')
  plt.plot([25.,25.],[ymin,ymax],c='grey',ls='dotted')
  plt.plot([0.,75.],[0.5*maxval,0.5*maxval],c='grey',ls='dotted')
  plt.legend(title='Aqueous',loc='upper right',fontsize=legend_fontsize)
  legend = plt.gca().get_legend()
  legend.get_frame().set_fill(False)
  legend.draw_frame(False)

aq_labels = []
aq_labels.append('Aaq')
aq_labels.append('Baq')
aq_labels.append('Caq')
aq_labels.append('Daq')
aq_labels.append('Eaq')
aq_labels.append('Faq')

aq_colors = []
aq_colors.append('blue')
aq_colors.append('green')
aq_colors.append('red')
aq_colors.append('cyan')
aq_colors.append('magenta')
aq_colors.append('y')

aq_linestyles = []
aq_linestyles.append('-')
aq_linestyles.append('-')
aq_linestyles.append('-')
aq_linestyles.append('--')
aq_linestyles.append('-')
aq_linestyles.append('--')

im_labels = []
im_labels.append('Xim')
im_labels.append('Yim')

im_colors = []
im_colors.append('darkorange')
im_colors.append('navy')

f = plt.figure(figsize=(16,6))
f.suptitle("Solute Transport, Sorption and Radioactive Decay",fontsize=16)

#linear scale
plt.subplot(1,2,1)
plt.title('Concentration Breakthrough @ 49.5 Meters')
scale_string = 'linear'
plot_aq_conc_breakthrough(plt,observation_filename,scale_string)

#log scale
plt.subplot(1,2,2)
plt.title('Concentration Profile @ 12.5 Years')
scale_string = 'linear'
plot_aq_conc_profile(plt,plot_filename,scale_string)

f.subplots_adjust(hspace=0.2,wspace=0.40,
                  bottom=.12,top=.85,
                  left=.08,right=.92)

plt.show()
