#Turn a Vorocrust .vcg mesh into a pflotran EXPLICIT_UNSTRUCTURED mesh
# and associated material ID and boundary files.
#Author: Tara LaForce 7/7/2021

#SYNTAX:
# python voro2pflotran.py mesh.vcg Xmin Xmax Ymin Ymax Zmin Zmax
#where Xmin-Zmax are the range of the meshed domain

#This script converts a vorocrust mesh.vcg file into a series of PFLOTRAN-readable files:
#1) out_mesh.uge is a PFLOTRAN explicit unstructured mesh.
#2) Seven boundary condition files are created:
#  a) bdry_smX.ex, bdry_lgX.ex, bdry_smY.ex, bdry_lgY.ex, bdry_smZ.ex, and
# bdry_lgZ.ex contain the boundary faces for the user-defined small and large x,y, and z
# boundaries if there are model boundaries aligned with the coordinate axes. If the model
# boundaries are not aligned with the axes or the incorrect x, y and z  boundaries are
# input these files are generated but contain no connections.
#  b) bdry_unclaim.ex contains all boundary elements that are not assigned in the previous
# step.  This file may also be used as a boundary condition if desired.
#3) matID.h5 contains the region material ID assigned by VoroCrust to each element.  It
# can be used to assign material IDs in a PFLOTRAN simulation.

#Examples:
#1) Use out_mesh.uge in PFLOTRAN as an explicit unstructured mesh
#GRID
#  TYPE unstructured_explicit out_mesh.uge
#  UPWIND_FRACTION_METHOD CELL_VOLUME
#END
#NOTE: It is highly recommended to use the UPWIND_FRACTION keyword for voronoi meshes

#2) Using bdry_*.ex file for boundary conditions.
#REGION innerX
#  FILE bdry_smX.ex
#END

#3) Read Material ID file matID.h5 into strata block.
#STRATA
#  FILE matID.h5
#END
#NOTE:A Material ID number is automatically assigned to every closed volume by VoroCrust
# and cannot be changed. Users need to change PFLOTRAN MATERIAL block to match the
# VoroCrust ID.

#NTESS statement:
#Sandia National Laboratories is a multimission laboratory managed and operated by
#National Technology & Engineering Solutions of Sandia, LLC, a wholly owned subsidiary
#of Honeywell International Inc., for the U.S. Department of Energy’s National Nuclear
#Security Administration under contract DE-NA0003525.

#Copyright 2021 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
#Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain
#rights in this software.

import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np
import sys
import h5py

def nan_image_connect(v1,v2):
    if (v1 <= n_cells) & (v2 <= n_cells):
        return v1
    else:
        return -99

print('**************************************************************************',\
      'Copyright 2021 National Technology and Engineering Solutions of Sandia, LLC.',\
      'Under the terms of Contract DE-NA0003525, there is a non-exclusive license',\
      'for use of this work by or on behalf of the U.S. Government. Export of this',\
      'program may require a license from the United States Government.',\
      '**************************************************************************',sep='\n')

if len(sys.argv) != 8:
    print("ERROR: Command line arguments not provided.\n")
    print("Usage: python load_voro_mesh_general.py input.vcg Xmin Xmax Ymin Ymax Zmin Zmax")
    sys.exit(0)
        
mesh_name = sys.argv[1]

bdry = np.zeros(6)
for i in range(0,6):
    bdry[i] = sys.argv[2+i]
print('Xmin Xmax Ymin Ymax Zmin Zmax for the boundary files')
print(bdry)

x_min = bdry[0]
x_max = bdry[1]
y_min = bdry[2]
y_max = bdry[3]
z_min = bdry[4]
z_max = bdry[5]


mesh_head=pd.read_csv(mesh_name,sep='\s+',header=None,index_col=None,nrows=1)
n_cells = mesh_head.iloc[0,1]
print('Number of cells in mesh')
print(n_cells)

mesh=pd.read_csv(mesh_name,sep='\s+',header=None,skiprows=1,index_col=0,nrows=n_cells)
mesh.columns = ['x','y','z','Volume','MatID']
#print(mesh.head())
#print(mesh.tail())
mesh_np=mesh.as_matrix()

file_name = 'out_'+ mesh_name.rstrip('.vcg') + '.uge'
print('New mesh file name is:')
print(file_name)

cell_num_str = ['CELLS '+ str(n_cells)]       
fd=open(file_name,'w')
np.savetxt(fd,cell_num_str,fmt='%1s')
for i in range(1,n_cells+1):
    this_row = np.append(i,mesh_np[i-1,0:4])
    np.savetxt(fd,this_row.reshape(1,this_row.shape[0]),fmt='%1i %1.15e %1.15e %1.15e %1.15e')
fd.close()


########################
connect=np.genfromtxt(mesh_name,skip_header= n_cells+2)

len_connect = len(connect[:,0])
real_connect = np.zeros((len_connect,6))
smX_connect = np.zeros((len_connect,5))
lgX_connect = np.zeros((len_connect,5))
smY_connect = np.zeros((len_connect,5))
lgY_connect = np.zeros((len_connect,5))
smZ_connect = np.zeros((len_connect,5))
lgZ_connect = np.zeros((len_connect,5))
unclaim_connect = np.zeros((len_connect,5))
reassigned_v1 = np.zeros((len_connect))
reassigned_v2 = np.zeros((len_connect))
num_real = 0
num_smX = 0
num_lgX = 0
num_smY = 0
num_lgY = 0
num_smZ = 0
num_lgZ = 0
num_unclaim = 0
del_x = 1e-8
print('number of connections read')
unclaimed_boundaries = 0
for i in range(0,len_connect):
    bdry_claimed = 0
    if math.ceil(i/10000)==math.floor(i/10000):
        print(i)
        print(num_real)
    if (connect[i,0] != connect[i,1]):
        real_connect[num_real,:]=connect[i,0:6]
        num_real = num_real + 1
        bdry_claimed = 1
        
    if (connect[i,2] > x_max-del_x) & (connect[i,2] < x_max+del_x) & (connect[i,0] == connect[i,1]):
        lgX_connect[num_lgX,0]=connect[i,0]
        lgX_connect[num_lgX,1:5]=connect[i,2:6]
        num_lgX = num_lgX + 1
        bdry_claimed = 1
    if (connect[i,2] > x_min-del_x) & (connect[i,2] < x_min+del_x) & (connect[i,0] == connect[i,1]):
        smX_connect[num_smX,0]=connect[i,0]
        smX_connect[num_smX,1:5]=connect[i,2:6]
        num_smX = num_smX + 1
        bdry_claimed = 1

    if (connect[i,3] > y_max-del_x) & (connect[i,3] < y_max+del_x) & (connect[i,0] == connect[i,1]):
        lgY_connect[num_lgY,0]=connect[i,0]
        lgY_connect[num_lgY,1:5]=connect[i,2:6]
        num_lgY = num_lgY + 1
        bdry_claimed = 1
    if (connect[i,3] > y_min-del_x) & (connect[i,3] < y_min+del_x) & (connect[i,0] == connect[i,1]):
        smY_connect[num_smY,0]=connect[i,0]
        smY_connect[num_smY,1:5]=connect[i,2:6]
        num_smY = num_smY + 1
        bdry_claimed = 1


    if (connect[i,4] > z_max-del_x) & (connect[i,4] < z_max+del_x) & (connect[i,0] == connect[i,1]):
        #print 'large Z connection'
        #print connect[i,:]
        lgZ_connect[num_lgZ,0]=connect[i,0]
        lgZ_connect[num_lgZ,1:5]=connect[i,2:6]
        num_lgZ = num_lgZ + 1
        bdry_claimed = 1
    if (connect[i,4] > z_min-del_x) & (connect[i,4] < z_min+del_x) & (connect[i,0] == connect[i,1]):
        smZ_connect[num_smZ,0]=connect[i,0]
        smZ_connect[num_smZ,1:5]=connect[i,2:6]
        num_smZ = num_smZ + 1
        bdry_claimed = 1
    if bdry_claimed == 0:
        print('ERROR! Boundary element has not been assigned to any boundary')
        print(connect[i,:])
        unclaim_connect[num_unclaim,0]=connect[i,0]
        unclaim_connect[num_unclaim,1:5]=connect[i,2:6]
        num_unclaim = num_unclaim + 1
        unclaimed_boundaries = unclaimed_boundaries + 1

print('Number of unclaimed boundaries')
print(unclaimed_boundaries)
fd=open(file_name,'a')
conne = ['CONNECTIONS '+ str(num_real)]       
np.savetxt(fd,conne,fmt='%5s')
np.savetxt(fd,real_connect[0:num_real,:],fmt='%1i %1i %1.15e %1.15e %1.15e %1.15e')
fd.close()

fd=open('bdry_smX.ex','w')
conne = ['CONNECTIONS '+ str(num_smX)]       
np.savetxt(fd,conne,fmt='%5s')
np.savetxt(fd,smX_connect[0:num_smX,:],fmt='%1i %1.15e %1.15e %1.15e %1.15e')
fd.close()


fd=open('bdry_lgX.ex','w')
conne = ['CONNECTIONS '+ str(num_lgX)]       
np.savetxt(fd,conne,fmt='%5s')
np.savetxt(fd,lgX_connect[0:num_lgX,:],fmt='%1i %1.15e %1.15e %1.15e %1.15e')
fd.close()

fd=open('bdry_smY.ex','w')
conne = ['CONNECTIONS '+ str(num_smY)]       
np.savetxt(fd,conne,fmt='%5s')
np.savetxt(fd,smY_connect[0:num_smY,:],fmt='%1i %1.15e %1.15e %1.15e %1.15e')
fd.close()


fd=open('bdry_lgY.ex','w')
conne = ['CONNECTIONS '+ str(num_lgY)]       
np.savetxt(fd,conne,fmt='%5s')
np.savetxt(fd,lgY_connect[0:num_lgY,:],fmt='%1i %1.15e %1.15e %1.15e %1.15e')
fd.close()

fd=open('bdry_smZ.ex','w')
conne = ['CONNECTIONS '+ str(num_smZ)]       
np.savetxt(fd,conne,fmt='%5s')
np.savetxt(fd,smZ_connect[0:num_smZ,:],fmt='%1i %1.15e %1.15e %1.15e %1.15e')
fd.close()


fd=open('bdry_lgZ.ex','w')
conne = ['CONNECTIONS '+ str(num_lgZ)]       
np.savetxt(fd,conne,fmt='%5s')
np.savetxt(fd,lgZ_connect[0:num_lgZ,:],fmt='%1i %1.15e %1.15e %1.15e %1.15e')
fd.close()

fd=open('bdry_unclaim.ex','w')
conne = ['CONNECTIONS '+ str(num_unclaim)]       
np.savetxt(fd,conne,fmt='%5s')
np.savetxt(fd,unclaim_connect[0:num_unclaim,:],fmt='%1i %1.15e %1.15e %1.15e %1.15e')
fd.close()

h5out = h5py.File('matID.h5', 'w')
#dset = f.create_dataset("mydataset", (100,), dtype='i')
len_mesh = n_cells
cell_ids = np.arange(1,len_mesh+1,dtype=int)

#por = np.array(h5['Checkpoint/PMCSubsurfaceFlow/flow/Porosity'][:],'=f8')
#h5['Checkpoint/PMCSubsurfaceFlow/flow/Porosity'][:] = por
  
h5out.create_dataset("Materials/Material Ids", data=mesh_np[:,4])
#mesh.matID)
h5out.create_dataset("Materials/Cell Ids", data=cell_ids)

h5out.close()
