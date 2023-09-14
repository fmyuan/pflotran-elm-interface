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

data_t1_out = pft.Dataset(filename,1,18)
data_t1_in = pft.Dataset(filename,1,5)
data_t2_out = pft.Dataset(filename,1,19)
data_t2_in = pft.Dataset(filename,1,6)
data_t3_out = pft.Dataset(filename,1,20)
data_t3_in = pft.Dataset(filename,1,7)

plt.semilogx(data_t1_out.get_array('x'),-data_t1_out.get_array('y')/data_t1_in.get_array('y')[0])
plt.semilogx(data_t2_out.get_array('x'),-data_t2_out.get_array('y')/data_t2_in.get_array('y')[0])
plt.semilogx(data_t3_out.get_array('x'),-data_t3_out.get_array('y')/data_t3_in.get_array('y')[0])
plt.xlabel("Time [yrs]")
plt.ylabel("Normalized Tracer")

plt.show()
