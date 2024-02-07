import sys
import h5py
sys.path.append(sys.path[0]+"/../../../../src/python/unstructured_grid")

from rotate_ugrid_implicit import *
from struct_grid_to_ugrid_implicit import *

supp = False

#generate reference grid
ref_name = "mesh_ref_ugi.h5"
h5 = True
n = (9,9,8) #nx,ny,nz
d = [[100.],[100.],[50.]] #uniform grid spacing in all direction
origin = (0,0,0)
struct_grid_to_ugrid_implicit(ref_name, h5, n, d, origin, True)


#def rotation to be done
rot1 = (23.,0., 0) #(angleXY, angleXZ, angleYZ)
rot2 = (0.,108.,0.)
rot3 = (0.,0.,-54.)
rot4 = (10.,36.,-75.)
rot_origin = (0.,0.,0.)

#rotate the grid
if supp:
  #rot1
  rot1_name = "mesh_rot1_ugi.h5"
  rotate_ugrid(ref_name, rot1_name, rot1, rot_origin)
  #rot2
  rot2_name = "mesh_rot2_ugi.h5"
  rotate_ugrid(ref_name, rot2_name, rot2, rot_origin)
  #rot3
  rot3_name = "mesh_rot3_ugi.h5"
  rotate_ugrid(ref_name, rot3_name, rot3, rot_origin)
#rot4
rot4_name = "mesh_rot4_ugi.h5"
rotate_ugrid(ref_name, rot4_name, (rot4[0],0.,0.), rot_origin) #xy first
rotate_ugrid(rot4_name, rot4_name, (0.,rot4[1],0.), rot_origin) #xz second
rotate_ugrid(rot4_name, rot4_name, (0.,0.,rot4[2]), rot_origin) #yz last


#create the groups
if supp:
  l = (ref_name, rot1_name, rot2_name, rot3_name, rot4_name)
else:
  l = (ref_name, rot4_name)

for x in l:
  f = h5py.File(x, 'r+')
  f.create_dataset('Regions/injection/Cell Ids', data=np.array([50]))
  f.create_dataset('Regions/extraction/Cell Ids', data=np.array([600]))
  f.create_dataset('Regions/obs1/Cell Ids', data=np.array([200]))
  f.create_dataset('Regions/obs2/Cell Ids', data=np.array([300]))
  f.create_dataset('Regions/obs3/Cell Ids', data=np.array([500]))
  f.close()
  
print_tensor_comp = False
if print_tensor_comp:
  #compute the tensor component
  kxx = 1e-13; kyy = 5e-14; kzz = 1e-14
  K_ref = (kxx,0,0,kyy,0,kzz)
  if supp:
    #rot1
    K1 = compute_perm_tensor_component(K_ref, rot1)
    print(K1)
    #rot2
    K2 = compute_perm_tensor_component(K_ref, rot2)
    print(K2)
    #rot3
    K3 = compute_perm_tensor_component(K_ref, rot3)
    print(K3)
  #rot4
  K4 = compute_perm_tensor_component(K_ref, (rot4[0],0.,0.))
  K4 = compute_perm_tensor_component(K4, (0.,rot4[1],0.))
  K4 = compute_perm_tensor_component(K4, (0.,0.,rot4[2]))
  print(K4)
