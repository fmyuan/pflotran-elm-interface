import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
import pflotran as pft

legend_fontsize = 'small'

default_labels = []
default_labels.append('Aaq (conservative)')
default_labels.append('Baq (sorption)')
default_labels.append('Caq (decay reactant)')
default_labels.append('Daq (decay product)')
default_labels.append('Eaq (sorption and decay reactant)')
default_labels.append('Faq (decay product)')

default_colors = []
default_colors.append('blue')
default_colors.append('green')
default_colors.append('red')
default_colors.append('cyan')
default_colors.append('magenta')
default_colors.append('y')

default_linestyles = ['-']*len(default_colors)

# legend location
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

def plot_results(plt,filenames,columns,
                 labels=default_labels,
                 colors=default_colors,
                 linestyles=default_linestyles,
                 scale_string='linear'):
    maxval = -1.e20
    minval = 1.e-10
    i = 0
    for filename in filenames:
        for column in columns:
            data = pft.Dataset(filename,1,column)
            ydata = data.get_array('y')
            maxval = max(maxval,np.amax(ydata))
            plt.plot(data.get_array('x'),data.get_array('y'),
                     label=labels[i],c=colors[i],ls=linestyles[i])
            i += 1
    if scale_string == 'linear':
        plt.ylim(-0.02*maxval,1.05*maxval)
    else:
        plt.ylim(0.5*minval,2.*maxval)
    plt.legend(loc='best',fontsize=legend_fontsize)
    legend = plt.gca().get_legend()
    legend.get_frame().set_fill(False)
    legend.draw_frame(False)
    return minval, maxval

def gridlines(plt,x,y):
    plt.plot(x,y,c='grey',ls='dotted')

def create_figure(title):
    f = plt.figure(figsize=(12,6))
    f.suptitle(title,fontsize=16)
    f.subplots_adjust(hspace=0.2,wspace=0.25,
                      bottom=.12,top=.85,
                      left=.1,right=.98)
    return f
