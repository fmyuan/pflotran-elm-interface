#
# Convert a FEHM mesh to PFLOTRAN unstructured implicit hdf5 format
#

import sys
import h5py
import numpy as np

def convert_fehm_grid_to_h5(fehm_in, h5_out):
    src = open(fehm_in,'r')
    #skip 'coor' keyword
    line = src.readline()
    num_vert = int(src.readline().strip())
    vertex_coordinates = np.zeros((num_vert,3),dtype='f8')
    for ivert in range(num_vert):
        line = src.readline()
        w = line.split()
        vertex_id = int(w[0])-1
        for i in range(3):
            vertex_coordinates[vertex_id][i] = float(w[i+1])
    # skip last line with 'int'
    line = src.readline()

    #skip 'elem' keyword
    line = src.readline()
    w = src.readline().split()
    num_elem = int(w[1])
    num_vertex_per_cell = int(w[0])
    element_vertices = np.zeros((num_elem,num_vertex_per_cell+1),dtype='i4')
    for ielem in range(num_elem):
        line = src.readline()
        w = line.split()
        cell_id = int(w[0])-1
        element_vertices[cell_id][0] = num_vertex_per_cell
        for i in range(1,num_vertex_per_cell+1):
            element_vertices[cell_id][i] = int(w[i])
        i = element_vertices[cell_id][2]
        element_vertices[cell_id][2] = element_vertices[cell_id][3]
        element_vertices[cell_id][3] = i
    src.close()
  
    #write to hdf5 output
    out = h5py.File(h5_out,'w')
    out.create_dataset("Domain/Cells", data=element_vertices)
    out.create_dataset("Domain/Vertices", data=vertex_coordinates)
    out.close()
    return 0


if __name__ == "__main__":
  if len(sys.argv) != 3:
    print("Utilization:")
    print("{} [path to FEHM grid] [output file]\n".format(sys.argv[0]))
    exit(1)
  fehm_in = sys.argv[1]
  h5_out = sys.argv[2]
  ret_code = convert_fehm_grid_to_h5(fehm_in, h5_out)
  exit(ret_code)
