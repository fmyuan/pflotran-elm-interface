import numpy as np
import pandas as pd
import re
import math
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import warnings

warnings.simplefilter('ignore')
# Do some crude post-processing

injector_file = 'injector.well'
producer_file = 'producer.well'

inj_copy = open(injector_file,'r')
start_string = 'Time [y]'
for i,line in enumerate(inj_copy):
    if start_string.lower() in line.lower():
        skiprows = i
        break
inj_copy.close()
with open(injector_file) as inj:
    for i in range(skiprows):
        inj.readline()
    while True:
        row = inj.readline().strip()
        if row:
            headers = re.split(r',',row)
            break
i = 0
num_inj_segs = 0
num_repeated_vars = 17
new_headers = [headers[0][1:-1],headers[1][1:-1]]
for word in headers:
    if i > 1:
        seg_number = math.ceil((i-1)/num_repeated_vars)
        word=word[1:-1]
        word+='_'+str(seg_number)
        new_headers.append(word)
        num_inj_segs += 1
    i += 1
headers = new_headers
injector = pd.read_csv(injector_file,skiprows=skiprows+1,sep='   ',header=None,names=headers)

prod_copy = open(producer_file,'r')
start_string = 'Time [y]'
for i,line in enumerate(prod_copy):
    if start_string.lower() in line.lower():
        skiprows = i
        break
prod_copy.close()
with open(producer_file) as inj:
    for i in range(skiprows):
        inj.readline()
    while True:
        row = inj.readline().strip()
        if row:
            headers = re.split(r',',row)
            break
i = 0
num_prod_segs = 0
num_repeated_vars = 17
new_headers = [headers[0][1:-1],headers[1][1:-1]]
for word in headers:
    if i > 1:
        seg_number = math.ceil((i-1)/num_repeated_vars)
        word=word[1:-1]
        word+='_'+str(seg_number)
        new_headers.append(word)
        num_prod_segs += 1
    i += 1

num_inj_segs = int(num_inj_segs/num_repeated_vars)
num_prod_segs = int(num_prod_segs/num_repeated_vars)
headers = new_headers
producer = pd.read_csv(producer_file,skiprows=skiprows+1,sep='   ',header=None,names=headers)

# Make plots
sec2yr = 3600*24*365
injection_start = 0.0
injection_end = 30.0
production_start = 0.0
production_end = 30.0
fig,(ax1,ax2) = plt.subplots(2,1)
fig.set_figheight(10)
fig.set_figwidth(10)
injector_filtered_gas = injector.filter(like='Well Q-gas [kg/s]').sum(axis=1)
injector_filtered_liq = injector.filter(like='Well Q-liq [kg/s]').sum(axis=1)
producer_filtered_gas = producer.filter(like='Well Q-gas [kg/s]').sum(axis=1)
producer_filtered_liq = producer.filter(like='Well Q-liq [kg/s]').sum(axis=1)
ax1.plot(injector['Time [y]'],-1*sec2yr * (injector_filtered_gas + injector_filtered_liq),'g--')
ax1.plot(producer['Time [y]'],-1*sec2yr * (producer_filtered_gas + producer_filtered_liq),'c--')
ax1.set_ylabel('Well Flow Rate [kg/y]')
ax1.set_xscale('linear')
ax1.set_xlim((1, 30))
ax1.spines['left'].set_visible(False)
ax1.yaxis.set_ticks_position('left')
ax1.yaxis.set_visible(False)
divider = make_axes_locatable(ax1)
ax1log=divider.append_axes("left", size=2.0, pad=0, sharey=ax1)
ax1log.set_xscale('log')
ax1log.set_xlim((1e-3, 1))
ax1log.plot(injector['Time [y]'],-1*sec2yr * (injector_filtered_gas + injector_filtered_liq),'g--')
ax1log.plot(producer['Time [y]'],-1*sec2yr * (producer_filtered_gas + producer_filtered_liq),'c--')
ax1log.spines['right'].set_visible(False)
plt.setp(ax1.get_xticklabels(), visible=True)

injector = injector[injector['Time [y]']<=injection_end][injector['Time [y]']>injection_start]
producer = producer[producer['Time [y]']<=production_end][producer['Time [y]']>production_start]
ax2.plot(injector['Time [y]'],injector['Well BHP [Pa]']/1e6)
ax2.plot(producer['Time [y]'],producer['Well BHP [Pa]']/1e6)
ax2.set_ylabel('Well Flow Rate [kg/y]')
ax2.set_xscale('linear')
ax2.set_xlim((1, 30))
ax2.spines['left'].set_visible(False)
ax2.yaxis.set_ticks_position('left')
ax2.yaxis.set_visible(False)
divider = make_axes_locatable(ax2)
ax2log=divider.append_axes("left", size=2.0, pad=0, sharey=ax2)
ax2log.set_xscale('log')
ax2log.set_xlim((1e-3, 1))
ax2log.plot(injector['Time [y]'],injector['Well BHP [Pa]']/1e6)
ax2log.plot(producer['Time [y]'],producer['Well BHP [Pa]']/1e6)
ax2log.spines['right'].set_visible(False)
plt.setp(ax2.get_xticklabels(), visible=True)

ax1.set_xlabel('Time (yr)')
ax1log.set_ylabel('Well Flow Rate (kg/yr)')
ax1.legend(['Injection Well','Production Well'])
ax2.set_xlabel('Time (yr)')
ax2log.set_ylabel('Well Pressure (MPa)')
ax2.legend(['Injection Well','Production Well'])
plt.savefig('multi-well.png')
