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
from plot_aux import *

f = create_figure("Solute Transport and Decay")

plt.subplot(1,2,1)
plt.title('Concentration Breakthrough @ 49.5 Meters')
plt.xlabel('Time [y]')
plt.ylabel('Concentration [M]')
filenames = ['transport_decay-obs-0.pft']
columns = [2,3,4]
mylabels = default_labels[:]
mylabels[1:2] = default_labels[2:4]
mycolors = default_colors[:]
mycolors[1:2] = default_colors[2:4]
minval,maxval = plot_results(plt,filenames,columns,
                             colors=mycolors,labels=mylabels)


plt.subplot(1,2,2)
plt.title('Concentration Profile @ 12.5 Years')
plt.xlabel('X [m]')
plt.ylabel('Concentration [M]')
filenames = ['transport_decay-002.tec']
columns = [4,5,6]
minval,maxval = plot_results(plt,filenames,columns,
                             colors=mycolors,labels=mylabels)

plt.show()
