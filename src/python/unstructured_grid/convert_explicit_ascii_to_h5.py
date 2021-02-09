#
# Convert an explicit unstructured grid in ascii to hdf5 format
#

import sys
import h5py
import numpy as np

def convert_grid_explicit_ascii_to_h5(ascii_in, h5_out):
  src = open(ascii_in,'r')
  #skip comment
  line = src.readline()
  while line[0] == '#': line = src.readline()
  #number of cells
  n_cells = int(line.split()[1])
  #initialize array
  cell_centers = np.zeros((n_cells,3),dtype='f8')
  cell_volumes = np.zeros(n_cells,dtype='f8')
  #read cell attributes (center and volume)
  for i in range(n_cells):
    line = [float(x) for x in src.readline().split()[1:]]
    cell_centers[i] = line[:3]
    cell_volumes[i] = line[3]
  #number of connections
  line = src.readline()
  n_connections = int(line.split()[1])
  #initialize array
  face_connections = np.zeros((n_connections,2),dtype='i8')
  face_areas = np.zeros(n_connections,dtype='f8')
  face_centers = np.zeros((n_connections,3),dtype='f8')
  # read face attributes (center, connections, area)
  for i in range(n_connections):
    line = src.readline().split()
    face_connections[i] = [int(x) for x in line[:2]]
    face_centers[i] = [float(x) for x in line[2:5]]
    face_areas[i] = float(line[5])
  src.close()
  
  #write to hdf5 output
  out = h5py.File(h5_out,'w')
  out.create_dataset("Domain/Cells/Centers", data=cell_centers)
  out.create_dataset("Domain/Cells/Volumes", data=cell_volumes)
  out.create_dataset("Domain/Connections/Centers", data=face_centers)
  out.create_dataset("Domain/Connections/Cell Ids", data=face_connections)
  out.create_dataset("Domain/Connections/Areas", data=face_areas)
  out.close()
  return 0


if __name__ == "__main__":
  if len(sys.argv) != 3:
    print("Utilisation:")
    print("{} [path to explicit ASCII grid] [output file]\n".format(sys.argv[0]))
    exit(1)
  ascii_in = sys.argv[1]
  h5_out = sys.argv[2]
  ret_code = convert_grid_explicit_ascii_to_h5(ascii_in, h5_out)
  exit(ret_code)
