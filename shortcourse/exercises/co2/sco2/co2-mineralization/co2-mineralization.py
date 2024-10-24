import numpy as np
import pandas as pd
import re
import math
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import warnings
from functools import reduce

warnings.simplefilter('ignore')

mass_balance_file = 'co2-mineralization-mas.dat'

with open(mass_balance_file) as mb:
    row = mb.readline().strip()
    headers = re.split(r',',row)
headers = reduce(lambda x, y: x + [re.sub('"', '', y)], headers, [])
mass_balance = pd.read_csv(mass_balance_file,skiprows=1,header=None,names=headers,sep='  | ')

# Make plots
f = plt.figure(figsize=(10,10))

plt.semilogy(mass_balance['Time [y]'],mass_balance['Region all Calcite Total Mass [mol]'],label='Calcite')
plt.semilogy(mass_balance['Time [y]'],mass_balance['Region all Anatase Total Mass [mol]'],label='Anatase')
plt.semilogy(mass_balance['Time [y]'],mass_balance['Region all Chalcedony Total Mass [mol]'],label='Chalcedony')
plt.semilogy(mass_balance['Time [y]'],mass_balance['Region all Dawsonite Total Mass [mol]'],label='Dawsonite')
plt.semilogy(mass_balance['Time [y]'],mass_balance['Region all Rhodochrosite Total M [mol]'],label='Rhodochrosite')

plt.xlabel('Time (yr)')
plt.ylabel('Mineral Mass (mol)')
plt.legend()

plt.savefig('minerals-formed.png')

f = plt.figure(figsize=(10,10))
plt.plot(mass_balance['Time [y]'],mass_balance['Region all Forsterite Total Mass [mol]']/mass_balance['Region all Forsterite Total Mass [mol]'][0],label='Forsterite')
plt.plot(mass_balance['Time [y]'],mass_balance['Region all Fayalite Total Mass [mol]']/mass_balance['Region all Fayalite Total Mass [mol]'][0],label='Fayalite')
plt.plot(mass_balance['Time [y]'],mass_balance['Region all Tephroite Total Mass [mol]']/mass_balance['Region all Tephroite Total Mass [mol]'][0],label='Tephroite')
plt.plot(mass_balance['Time [y]'],mass_balance['Region all Wollastonite Total Ma [mol]']/mass_balance['Region all Wollastonite Total Ma [mol]'][0],label='Wollastonite')
plt.plot(mass_balance['Time [y]'],mass_balance['Region all Enstatite Total Mass [mol]']/mass_balance['Region all Enstatite Total Mass [mol]'][0],label='Enstatite')
plt.plot(mass_balance['Time [y]'],mass_balance['Region all Chrysotile Total Mass [mol]']/mass_balance['Region all Chrysotile Total Mass [mol]'][0],label='Chrysotile')

plt.xlabel('Time (yr)')
plt.ylabel('Fraction of Initial Mineral Mass $\mathregular{M(t)/M_0}$')
plt.legend()

plt.savefig('mineral-change.png')

f = plt.figure(figsize=(10,10))
plt.plot(mass_balance['Time [y]'],mass_balance['Global CO2 Mass in Water Phase [kg]']/1e9,label='Aqueous')
plt.plot(mass_balance['Time [y]'],mass_balance['Global CO2 Mass in Gas Phase [kg]']/1e9,label='Free-Phase')
plt.plot(mass_balance['Time [y]'],mass_balance['Well well1 Total Gas Mass [kg]']/1e9,label='Total Injected')
plt.plot(mass_balance['Time [y]'],(mass_balance['Well well1 Total Gas Mass [kg]'] - mass_balance['Global CO2 Mass in Gas Phase [kg]'] - mass_balance['Global CO2 Mass in Water Phase [kg]'])/1e9,label='Mineralized')

plt.xlabel('Time (yr)')
plt.ylabel('CO2 Mass (MMT)')
plt.legend()
plt.savefig('co2-mass.png')
