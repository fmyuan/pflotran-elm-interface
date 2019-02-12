import sys
import os
try:
    pflotran_dir = os.environ['PFLOTRAN_DIR']
except KeyError:
    print('PFLOTRAN_DIR must point to PFLOTRAN installation directory and be defined in system environment variables.')
    sys.exit(1)
sys.path.append(pflotran_dir + '/src/python')

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
import pflotran as pft

mpl.rcParams['font.size']=12
#mpl.rcParams['font.weight']='bold'
mpl.rcParams['lines.linewidth']=2

def make_subplots(nrow, ncol, nplots, filenames, columns, name, titles):
  for iplot in range(nplots):
    plt.subplot(nrow,ncol,iplot+1)
    plt.title(titles[iplot])
    plt.xlabel('Time (years)', fontsize=14, fontweight='bold')
    plt.ylabel('I-129 Concentration (mol/L)', fontsize=14, fontweight='bold')

    plt.xlim(0.,5.e5)
    plt.ylim(1.e-22,1.e-9)
    plt.grid(True)

    for icol in range(iplot*len(columns)/nplots,iplot*len(columns)/nplots+len(columns)/nplots):
      ifile = icol
      data = pft.Dataset(filenames[ifile],1,columns[icol])
      string = name[icol]
      plt.semilogy(data.get_array('x'),data.get_array('y'),color=colors[icol],linestyle=linestyles[icol],label=string)

path = []
path.append('.')

titles = []
files = []
columns = []
name = []
linestyles = []
colors = []
basestyles = ['-','--'] #'-.',':'
basecolors = ['red','blue']

#Upper aquifer obs points
titles.append('Shallow Sandstone Aquifer')
files.append('./gdsa-obs-0.tec')
columns.append(29)
name.append('Near')
files.append('./gdsa-obs-0.tec')
columns.append(41)
name.append('Far')

linestyles.extend(basestyles)
colors.extend(basecolors)

#Shale host rock obs points
titles.append('Shale Host Rock')
files.append('./gdsa-obs-0.tec')
columns.append(5)
name.append('Near')
files.append('./gdsa-obs-0.tec')
columns.append(17)
name.append('Far')

linestyles.extend(basestyles)
colors.extend(basecolors)

#Lower aquifer obs points
titles.append('Deep Limestone Aquifer')
files.append('./gdsa-obs-0.tec')
columns.append(53)
name.append('Near')
files.append('./gdsa-obs-0.tec')
columns.append(65)
name.append('Far')

linestyles.extend(basestyles)
colors.extend(basecolors)

filenames = pft.get_full_paths(path,files)

f = plt.figure(figsize=(24,8)) #adjust size to fit plots (width, height)

make_subplots(1,3,3,filenames,columns,name,titles)

#Choose a location for the legend. #This is associated with the last plot
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
#plt.legend(loc=(1.20,0.175),title='FRACTURE REALIZATION')
plt.legend(loc=0)
# xx-small, x-small, small, medium, large, x-large, xx-large, 12, 14
plt.setp(plt.gca().get_legend().get_texts(),fontsize=14,fontweight='normal')
plt.setp(plt.gca().get_legend().get_texts(),linespacing=0.)
plt.setp(plt.gca().get_legend().get_frame().set_fill(False))
plt.setp(plt.gca().get_legend().draw_frame(False))
#plt.gca().yaxis.get_major_formatter().set_powerlimits((-1,1))

#adjust blank space and show/save the plot
f.subplots_adjust(hspace=0.2,wspace=0.2, bottom=.12,top=.9, left=.14,right=.9)
plt.savefig('gdsa_i129.png',bbox_inches='tight')
plt.show()
