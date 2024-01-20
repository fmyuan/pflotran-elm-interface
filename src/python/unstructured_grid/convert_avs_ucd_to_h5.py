#
# Convert an AVS UCD mesh to PFLOTRAN unstructured implicit hdf5 format
#

import sys
import h5py
import numpy as np

map_tet = [1,2,4,3]
map_pyr = [5,1,2,3,4]
map_prism = [4,5,6,1,2,3]
map_hex = [5,6,7,8,1,2,3,4]

def convert_avs_ucd_to_h5(ucd_in, h5_out):
    src = open(ucd_in,'r')

    w = src.readline().strip().split()
    num_vert = int(w[0])
    num_cell = int(w[1])
    vertex_coordinates = np.zeros((num_vert,3),dtype='f8')
    increment = max(int(num_vert/10),1)
    for ivert in range(num_vert):
        line = src.readline()
        w = line.split()
        vertex_id = int(w[0])-1
        for i in range(3):
            vertex_coordinates[vertex_id][i] = float(w[i+1])
        if ivert % increment == 0:
            print(' {:7d} vertices read'.format(ivert+1))
    print(' {:7d} vertices read'.format(num_vert))

    cell_vertices_pft = np.zeros((num_cell,9),dtype='i4')
    increment = max(int(num_cell/10),1)
    xmf_vertex_count = 0
    for icell in range(num_cell):
        line = src.readline()
        w = line.split()
        cell_id = int(w[0])-1
        cell_type = w[2]
        if cell_type == 'hex':
            map_ = map_hex
        elif cell_type == 'tet':
            map_ = map_tet
        elif cell_type == 'pyr':
            map_ = map_pyr
        elif cell_type == 'prism':
            map_ = map_prism
        else:
          print('Unrecognized AVS UCD cell type: {}'.format(cell_type))
          sys.exit(0)
        cell_vertices_pft[cell_id][0] = len(map_)
        for i in range(cell_vertices_pft[cell_id][0]):
            # map already includes +1 offset
            cell_vertices_pft[cell_id][map_[i]] = int(w[i+3])
        xmf_vertex_count += len(map_)
        if icell % increment == 0:
            print(' {:7d} cells read'.format(icell+1))
    print(' {:7d} cells read'.format(num_cell))
    src.close()
  
    #write to pflotran hdf5 output
    h5_prefix = h5_out.strip().rsplit('.')[0]
    out = h5py.File(h5_prefix+'_ugi.h5','w')
    out.create_dataset("Domain/Cells", data=cell_vertices_pft)
    out.create_dataset("Domain/Vertices", data=vertex_coordinates)
    out.close()
    #write to xmf hdf5 output
    out = h5py.File(h5_prefix+'_xmf.h5','w')
    cell_vertices_xmf = np.zeros((num_cell+xmf_vertex_count),dtype='i4')
    offset = 0
    for icell in range(num_cell):
        num_vert_in_cell = cell_vertices_pft[icell][0]
        if num_vert_in_cell == 8:
            xmf_cell_type = 9
        elif num_vert_in_cell == 4:
            xmf_cell_type = 6
        elif num_vert_in_cell == 5:
            xmf_cell_type = 7
        elif num_vert_in_cell == 6:
            xmf_cell_type = 8
        cell_vertices_xmf[offset] = xmf_cell_type
        offset += 1
        for ivert in range(num_vert_in_cell):
            # zero-based
            cell_vertices_xmf[offset] = cell_vertices_pft[icell][ivert+1]-1
            offset += 1
    if offset != num_cell+xmf_vertex_count:
        print('Incorrect number of entries in xmf Domain/Cells {} '.
              format(offset) + \
              ' when expecting {}'.format(num_cell+xmf_vertex_count))
        sys.exit(0)
    out.create_dataset("Domain/Cells", data=cell_vertices_xmf)
    out.create_dataset("Domain/Vertices", data=vertex_coordinates)
    out.close()

    return 0


if __name__ == "__main__":
  if len(sys.argv) != 3:
    print("Utilization:")
    print("{} [path to AVS UCS file] [output file]\n".format(sys.argv[0]))
    exit(1)
  ucd_in = sys.argv[1]
  h5_out = sys.argv[2]
  ret_code = convert_avs_ucd_to_h5(ucd_in, h5_out)
  exit(ret_code)
