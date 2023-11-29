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
sys.path.append(pflotran_dir + '/src/python/analytical')
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from ogata import *
from plot_aux import *

f = create_figure("Solute Transport")

plt.subplot(1,2,1)
plt.title('Concentration Breakthrough @ 49.5 Meters')
plt.xlabel('Time [y]')
plt.ylabel('Concentration [M]')
filenames = ['transport-obs-0.pft']
columns = [2]
minval,maxval = plot_results(plt,filenames,columns)

plt.subplot(1,2,2)
plt.title('Concentration Profile @ 12.5 Years')
plt.xlabel('X [m]')
plt.ylabel('Concentration [M]')
filenames = ['transport-002.tec']
columns = [4]
minval,maxval = plot_results(plt,filenames,columns)
gridlines(plt,[50.,50.],[minval,maxval])
gridlines(plt,[0.,100.],[0.5*maxval,0.5*maxval])

# ---------------------------------------------------------------------------- #
# analytical solutions
C0 = 1.e-3
velocity = 4. # m/y
coefficient_diffusion = 1.e-8
coefficient_diffusion *= 3600*24*365 # m^2/y
Deff = coefficient_diffusion
time = 12.5
distance = 49.5

# breakthrough
plt.subplot(1,2,1)
t = pft.Dataset('transport-obs-0.pft',1,2).get_array('x')
t[0] = 1.e-40 # avoid time = 0 that result in a divide by zero in ogata
nt = len(t)
c = np.zeros(nt)
for i in range(nt):
    c[i] = ogata(C0,Deff,velocity,distance,t[i])
plt.plot(t,c,ls='--',label='Analytical',c='black')
plt.legend(loc='best',fontsize=legend_fontsize)
legend = plt.gca().get_legend()
legend.get_frame().set_fill(False)
legend.draw_frame(False)

# concentration profile
plt.subplot(1,2,2)
x = pft.Dataset('transport-002.tec',1,2).get_array('x')
nx = len(x)
c = np.zeros(nx)
for i in range(nx):
    c[i] = ogata(C0,Deff,velocity,x[i],time)
plt.plot(x,c,ls='--',label='Analytical',c='black')
plt.legend(loc='best',fontsize=legend_fontsize)
legend = plt.gca().get_legend()
legend.get_frame().set_fill(False)
legend.draw_frame(False)

plt.show()
