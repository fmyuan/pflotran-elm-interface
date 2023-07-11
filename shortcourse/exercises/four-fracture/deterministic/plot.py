#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  4 15:19:26 2023

@author: rosie
"""

import sys,os
import pflotran as pft
import matplotlib.pyplot as plt

filename = "cpm_transport-mas.dat"

data_out = pft.Dataset(filename,1,12)
data_in = pft.Dataset(filename,1,5)

plt.plot(data_out.get_array('x'),-data_out.get_array('y')/data_in.get_array('y')[0])
plt.xlabel("Time [yrs]")
plt.ylabel("Normalized Tracer")
plt.show()
