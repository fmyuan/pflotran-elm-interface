import h5py
import numpy
from get_dimensions import get_dimensions

tet_map = numpy.zeros((6,4),'=i4')
tet_map[0,:] = [1,8,3,4]
tet_map[1,:] = [1,3,8,5]
tet_map[2,:] = [8,7,5,3]
tet_map[3,:] = [1,2,3,5]
tet_map[4,:] = [2,5,7,3]
tet_map[5,:] = [2,7,5,6]

def struct_grid_to_ugrid_implicit(out_name, h5, n, d, origin, hex_):
    """
    Create a unstructured implicit grid for PFLOTRAN
    out_name: output file name
    h5: true for hdf5 output or false for ascii
    n = [nx, ny nz] with n_i number of element in the i direction
    d = [[dx],[dy],[dz]] with [d_i] the list spacing in the i direction 
                                  (for uniform, one value only for each d_i
    origin = (x_o,y_o,z_o) coordinate of the origin
    hex_: element type
    """
    if h5:
        out = h5py.File(out_name, mode = 'w')
    else:
        out = open(out_name, 'w')

    nx,ny,nz = n
    dx,dy,dz = d
    nxp1 = nx+1; nyp1 = ny+1; nzp1 = nz+1

    if hex_:
        element_array = numpy.zeros((nx*ny*nz,9),'=i4')
        icell = 0
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    element_array[icell,0] = 8
                    for ivert in range(1,9):
                        element_array[icell,ivert] = \
                            local_vertex_to_offset(i,j,k,nxp1,nyp1,ivert)
                    icell += 1
    else:
        element_array = numpy.zeros((nx*ny*nz*6,9),'=i4')
        icell = 0
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    for itet in range(6):
                        element_array[icell,0] = 4
                        for ii in range(4):
                            element_array[icell,ii+1] = \
                                local_vertex_to_offset(i,j,k,nxp1,nyp1,tet_map[itet,ii])
                        icell += 1

    vertex_array = numpy.zeros((nxp1*nyp1*nzp1,3),'=f8')
    x_origin,y_origin,z_origin = origin

    for k in range(nzp1):
        if len(dz) == 1:
            z = k*dz[0] + z_origin
        else:
            if k == 0:
                z = z_origin
            else:
                z += dz[k-1]
        for j in range(nyp1):
            if len(dy) == 1:
                y = j*dy[0] + y_origin
            else:
                if j == 0:
                    y = y_origin
                else:
                    y += dy[j-1]
            for i in range(nxp1):
                vertex_id = i + j*nxp1 + k*nxp1*nyp1
                if len(dx) == 1:
                    x = i*dx[0] + x_origin
                else:
                    if i == 0:
                        x = x_origin
                    else:
                        x += dx[i-1]
                vertex_array[vertex_id,0] = x
                vertex_array[vertex_id,1] = y
                vertex_array[vertex_id,2] = z

    if h5:
        h5dset = out.create_dataset('Domain/Cells', data=element_array)
        h5dset = out.create_dataset('Domain/Vertices', data=vertex_array)
    else:
        num_cells = nx*ny*nz
        if not hex_:
            num_cells *= 6
        out.write('%d %d\n' % (num_cells,nxp1*nyp1*nzp1))
        for icell in range(num_cells):
              if hex_:
                  out.write('H')
              else:
                  out.write('T')
              for iv in range(element_array[icell,0]):
                  out.write(' %d' % (element_array[icell,iv+1]+1))
              out.write('\n')
        for i in range(nxp1*nyp1*nzp1):
              out.write('%12.6e %12.6e %12.6e\n' %
                        (vertex_array[i,0],vertex_array[i,1],vertex_array[i,2]))

    out.close()

    if hex_:
        west_faces = []
        east_faces = []
        for k in range(nz):
            for j in range(ny):
                west_faces.append(cell_ijk_to_face_vertices(0,j,k,nxp1,nyp1,'west'))
                east_faces.append(cell_ijk_to_face_vertices(nx-1,j,k,nxp1,nyp1,'east'))
        south_faces = []
        north_faces = []
        for k in range(nz):
            for i in range(nx):
                south_faces.append(cell_ijk_to_face_vertices(i,0,k,nxp1,nyp1,'south'))
                north_faces.append(cell_ijk_to_face_vertices(i,ny-1,k,nxp1,nyp1,'north'))
        bottom_faces = []
        top_faces = []
        for j in range(ny):
            for i in range(nx):
                bottom_faces.append(cell_ijk_to_face_vertices(i,j,0,nxp1,nyp1,'bottom'))
                top_faces.append(cell_ijk_to_face_vertices(i,j,nz-1,nxp1,nyp1,'top'))

    if not h5:
        output_faces('mesh_ugi_west.ss',west_faces,hex_)
        output_faces('mesh_ugi_east.ss',east_faces,hex_)
        output_faces('mesh_ugi_south.ss',south_faces,hex_)
        output_faces('mesh_ugi_north.ss',north_faces,hex_)
        output_faces('mesh_ugi_bottom.ss',bottom_faces,hex_)
        output_faces('mesh_ugi_top.ss',top_faces,hex_)

    return True

def vertex_ijk_to_vertex_offset(i,j,k,nxp1,nyp1):
    # i,j,k are vertex indices
    return i + j*nxp1 + k*nxp1*nyp1

def cell_ijk_to_face_vertices(i,j,k,nxp1,nyp1,face):
    vertices = []
    if face == 'west':
        vertices.append(vertex_ijk_to_vertex_offset(i,j,k,nxp1,nyp1))
        vertices.append(vertex_ijk_to_vertex_offset(i,j,k+1,nxp1,nyp1))
        vertices.append(vertex_ijk_to_vertex_offset(i,j+1,k+1,nxp1,nyp1))
        vertices.append(vertex_ijk_to_vertex_offset(i,j+1,k,nxp1,nyp1))
    elif face == 'east':
        vertices.append(vertex_ijk_to_vertex_offset(i+1,j,k,nxp1,nyp1))
        vertices.append(vertex_ijk_to_vertex_offset(i+1,j+1,k,nxp1,nyp1))
        vertices.append(vertex_ijk_to_vertex_offset(i+1,j+1,k+1,nxp1,nyp1))
        vertices.append(vertex_ijk_to_vertex_offset(i+1,j,k+1,nxp1,nyp1))
    elif face == 'south':
        vertices.append(vertex_ijk_to_vertex_offset(i,j,k,nxp1,nyp1))
        vertices.append(vertex_ijk_to_vertex_offset(i+1,j,k,nxp1,nyp1))
        vertices.append(vertex_ijk_to_vertex_offset(i+1,j,k+1,nxp1,nyp1))
        vertices.append(vertex_ijk_to_vertex_offset(i,j,k+1,nxp1,nyp1))
    elif face == 'north':
        vertices.append(vertex_ijk_to_vertex_offset(i,j+1,k,nxp1,nyp1))
        vertices.append(vertex_ijk_to_vertex_offset(i,j+1,k+1,nxp1,nyp1))
        vertices.append(vertex_ijk_to_vertex_offset(i+1,j+1,k+1,nxp1,nyp1))
        vertices.append(vertex_ijk_to_vertex_offset(i+1,j+1,k,nxp1,nyp1))
    elif face == 'bottom':
        vertices.append(vertex_ijk_to_vertex_offset(i,j,k,nxp1,nyp1))
        vertices.append(vertex_ijk_to_vertex_offset(i,j+1,k,nxp1,nyp1))
        vertices.append(vertex_ijk_to_vertex_offset(i+1,j+1,k,nxp1,nyp1))
        vertices.append(vertex_ijk_to_vertex_offset(i+1,j,k,nxp1,nyp1))
    elif face == 'top':
        vertices.append(vertex_ijk_to_vertex_offset(i,j,k+1,nxp1,nyp1))
        vertices.append(vertex_ijk_to_vertex_offset(i+1,j,k+1,nxp1,nyp1))
        vertices.append(vertex_ijk_to_vertex_offset(i+1,j+1,k+1,nxp1,nyp1))
        vertices.append(vertex_ijk_to_vertex_offset(i,j+1,k+1,nxp1,nyp1))
    return vertices

def local_vertex_to_offset(i,j,k,nxp1,nyp1,vertex_id):
    # i,j,k are vertex indices
    offset = -999
    if vertex_id == 1:
        offset = vertex_ijk_to_vertex_offset(i,j,k,nxp1,nyp1)
    elif vertex_id == 2:
        offset = vertex_ijk_to_vertex_offset(i+1,j,k,nxp1,nyp1)
    elif vertex_id == 3:
        offset = vertex_ijk_to_vertex_offset(i+1,j+1,k,nxp1,nyp1)
    elif vertex_id == 4:
        offset = vertex_ijk_to_vertex_offset(i,j+1,k,nxp1,nyp1)
    elif vertex_id == 5:
        offset = vertex_ijk_to_vertex_offset(i,j,k+1,nxp1,nyp1)
    elif vertex_id == 6:
        offset = vertex_ijk_to_vertex_offset(i+1,j,k+1,nxp1,nyp1)
    elif vertex_id == 7:
        offset = vertex_ijk_to_vertex_offset(i+1,j+1,k+1,nxp1,nyp1)
    else:
        offset = vertex_ijk_to_vertex_offset(i,j+1,k+1,nxp1,nyp1)
    return offset

def output_faces(filename,faces,hex_):
    out = open(filename, 'w')
    out.write('{}\n'.format(len(faces)))
    for face_vertices in faces:
        if hex_:
            out.write('Q ')
        else:
            out.write('T ')
        for vertex in face_vertices:
            out.write('{} '.format(vertex+1))
        out.write('\n')
    out.close()

if __name__ == "__main__":

    h5 = False
    if h5:
        out_name = "mesh_ugi.h5"
    else:
        out_name = "mesh.ugi"

    hex_ = True

    n, d, origin = get_dimensions()
    struct_grid_to_ugrid_implicit(out_name, h5, n, d, origin, hex_)

    print('Done with unstructured implicit grid.')

