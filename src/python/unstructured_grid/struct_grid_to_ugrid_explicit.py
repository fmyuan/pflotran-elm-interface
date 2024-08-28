from h5py import *
import numpy
from get_dimensions import get_dimensions

def struct_grid_to_ugrid_explicit(out_name, h5, n, d, origin):
    """
    Create a unstructured explicit grid for PFLOTRAN.
    out_name: output file name
    h5: true for hdf5 output or false for ascii
    n = [nx, ny nz] with n_i number of element in the i direction
    d = [[dx],[dy],[dz]] with [d_i] the list spacing in the i direction 
                                  (for uniform, one value only for each d_i
    origin = (x_o,y_o,z_o) coordinate of the origin
    """
    f = open(out_name,'w')

    x_origin,y_origin,z_origin = origin
    nx,ny,nz = n
    dx,dy,dz = d

    x_face_offsets = []
    y_face_offsets = []
    z_face_offsets = []

    x = x_origin
    x_face_offsets.append(x)
    for i in range(nx):
        if len(dx) == 1:
            x += dx[0]
        else:
            x += dx[i]
        x_face_offsets.append(x)
    y = y_origin
    y_face_offsets.append(y)
    for j in range(ny):
        if len(dy) == 1:
            y += dy[0]
        else:
            y += dy[j]
        y_face_offsets.append(y)
    z = z_origin
    z_face_offsets.append(z)
    for k in range(nz):
        if len(dz) == 1:
            z += dz[0]
        else:
            z += dz[k]
        z_face_offsets.append(z)
    face_offsets = [x_face_offsets,y_face_offsets,z_face_offsets]

    f.write('CELLS %d\n' % (nx*ny*nz))
    z = z_origin + 0.5*dz[0]
    for k in range(nz):
        if len(dz) == 1:
            lenz = dz[0]
        else:
            lenz = dz[k]
        y = y_origin + 0.5*dy[0]
        for j in range(ny):
            if len(dy) == 1:
                leny = dy[0]
            else:
                leny = dy[j]
            x = x_origin + 0.5*dx[0]
            for i in range(nx):
                cell_id = i + j*nx + k*nx*ny
                if len(dx) == 1:
                    lenx = dx[0]
                else:
                    lenx = dx[i]
                volume = lenx*leny*lenz
                f.write('%d %f %f %f %f\n' % (cell_id + 1,x,y,z,volume))
                if len(dx) == 1:
                    x += dx[0]
                else:
                    x += 0.5*dx[i]
                    if i+1 < nx:
                        x += 0.5*dx[i+1]
            if len(dy) == 1:
                y += dy[0]
            else:
                y += 0.5*dy[j]
                if j+1 < ny:
                    y += 0.5*dy[j+1]
        if len(dz) == 1:
            z += dz[0]
        else:
            z += 0.5*dz[k]
            if k+1 < nz:
                z += 0.5*dz[k+1]

    f.write('CONNECTIONS %d\n' % ((nx-1)*ny*nz+nx*(ny-1)*nz+nx*ny*(nz-1)))

    # x connections
    z = z_origin + 0.5*dz[0]
    for k in range(nz):
        if len(dz) == 1:
            lenz = dz[0]
        else:
            lenz = dz[k]
        y = y_origin + 0.5*dy[0]
        for j in range(ny):
            if len(dy) == 1:
                leny = dy[0]
            else:
                leny = dy[j]
            x = x_origin + dx[0]
            for i in range(nx-1):
                cell_id = i + j*nx + k*nx*ny
                area = leny*lenz
                f.write('%d %d %f %f %f %f\n' % (cell_id + 1,
                                                 cell_id + 2,
                                                 x,y,z,area))
                if len(dx) == 1:
                    x += dx[0]
                else:
                    x += dx[i]
            if len(dy) == 1:
                y += dy[0]
            else:
                y += 0.5*dy[j]
                if j+1 < ny:
                    y += 0.5*dy[j+1]
        if len(dz) == 1:
            z += dz[0]
        else:
            z += 0.5*dz[k]
            if k+1 < nz:
                z += 0.5*dz[k+1]

    # y connections
    z = z_origin + 0.5*dz[0]
    for k in range(nz):
        if len(dz) == 1:
            lenz = dz[0]
        else:
            lenz = dz[k]
        x = x_origin + 0.5*dx[0]
        for i in range(nx):
            if len(dx) == 1:
                lenx = dx[0]
            else:
                lenx = dx[i]
            y = y_origin + dy[0]
            for j in range(ny-1):
                cell_id = i + j*nx + k*nx*ny
                area = lenx*lenz
                f.write('%d %d %f %f %f %f\n' % (cell_id + 1,
                                                 cell_id + 1 + nx,
                                                 x,y,z,area))
                if len(dy) == 1:
                    y += dy[0]
                else:
                    y += dy[j]
            if len(dx) == 1:
                x += dx[0]
            else:
                x += 0.5*dx[i]
                if i+1 < nx:
                    x += 0.5*dx[i+1]
        if len(dz) == 1:
            z += dz[0]
        else:
            z += 0.5*dz[k]
            if k+1 < nz:
                z += 0.5*dz[k+1]

    # z connections
    y = y_origin + 0.5*dy[0]
    for j in range(ny):
        if len(dy) == 1:
            leny = dy[0]
        else:
            leny = dy[j]
        x = x_origin + 0.5*dx[0]
        for i in range(nx):
            if len(dx) == 1:
                lenx = dx[0]
            else:
                lenx = dx[i]
            z = z_origin + dz[0]
            for k in range(nz-1):
                cell_id = i + j*nx + k*nx*ny
                area = leny*lenx
                f.write('%d %d %f %f %f %f\n' % (cell_id + 1,
                                                 cell_id + 1 + nx*ny,
                                                 x,y,z,area))
                if len(dz) == 1:
                    z += dz[0]
                else:
                    z += dz[k]
            if len(dx) == 1:
                x += dx[0]
            else:
                x += 0.5*dx[i]
                if i+1 < nx:
                    x += 0.5*dx[i+1]
        if len(dy) == 1:
            y += dy[0]
        else:
            y += 0.5*dy[j]
            if j+1 < ny:
                y += 0.5*dy[j+1]

    f.close()

    # boundaries
    west_faces = []
    east_faces = []
    for k in range(nz):
        for j in range(ny):
            west_faces.append(cell_ijk_to_face_connection(0,j,k,nx,ny,dx,dy,dz,
                                                          face_offsets,'west'))
            east_faces.append(cell_ijk_to_face_connection(nx-1,j,k,nx,ny,dx,dy,dz,
                                                          face_offsets,'east'))
    south_faces = []
    north_faces = []
    for k in range(nz):
        for i in range(nx):
            south_faces.append(cell_ijk_to_face_connection(i,0,k,nx,ny,dx,dy,dz,
                                                           face_offsets,'south'))
            north_faces.append(cell_ijk_to_face_connection(i,ny-1,k,nx,ny,dx,dy,dz,
                                                           face_offsets,'north'))
    bottom_faces = []
    top_faces = []
    for j in range(ny):
        for i in range(nx):
            bottom_faces.append(cell_ijk_to_face_connection(i,j,0,nx,ny,dx,dy,dz,
                                                            face_offsets,'bottom'))
            top_faces.append(cell_ijk_to_face_connection(i,j,nz-1,nx,ny,dx,dy,dz,
                                                         face_offsets,'top'))

    if not h5:
        output_face_connections('mesh_uge_west.ex',west_faces)
        output_face_connections('mesh_uge_east.ex',east_faces)
        output_face_connections('mesh_uge_south.ex',south_faces)
        output_face_connections('mesh_uge_north.ex',north_faces)
        output_face_connections('mesh_uge_bottom.ex',bottom_faces)
        output_face_connections('mesh_uge_top.ex',top_faces)

def cell_ijk_to_offset(i,j,k,nx,ny):
    return i + j*nx + k*nx*ny

def cell_ijk_to_face_connection(i,j,k,nx,ny,dxx,dyy,dzz,face_offsets,face):
    offset =  cell_ijk_to_offset(i,j,k,nx,ny)
    if len(dxx) == 1:
        dx = dxx[0]
    else:
        dx = dxx[i]
    if len(dyy) == 1:
        dy = dyy[0]
    else:
        dy = dyy[j]
    if len(dzz) == 1:
        dz = dzz[0]
    else:
        dz = dzz[k]
    if face == 'west':
        face_center_coord =[face_offsets[0][i],
                            face_offsets[1][j]+0.5*dy,
                            face_offsets[2][k]+0.5*dz]
        face_area = dy*dz
    elif face == 'east':
        face_center_coord =[face_offsets[0][i+1],
                            face_offsets[1][j]+0.5*dy,
                            face_offsets[2][k]+0.5*dz]
        face_area = dy*dz
    elif face == 'south':
        face_center_coord =[face_offsets[0][i]+0.5*dx,
                            face_offsets[1][j],
                            face_offsets[2][k]+0.5*dz]
        face_area = dx*dz
    elif face == 'north':
        face_center_coord =[face_offsets[0][i]+0.5*dx,
                            face_offsets[1][j+1],
                            face_offsets[2][k]+0.5*dz]
        face_area = dx*dz
    elif face == 'bottom':
        face_center_coord =[face_offsets[0][i]+0.5*dx,
                            face_offsets[1][j]+0.5*dy,
                            face_offsets[2][k]]
        face_area = dx*dy
    elif face == 'top':
        face_center_coord =[face_offsets[0][i]+0.5*dx,
                            face_offsets[1][j]+0.5*dy,
                            face_offsets[2][k+1]]
        face_area = dx*dy
    return [offset,face_center_coord,face_area]

def output_face_connections(filename,faces):
    out = open(filename, 'w')
    out.write('CONNECTIONS {}\n'.format(len(faces)))
    zero_to_one_based = 1
    for connection in faces:
        out.write('{} '.format(connection[0]+zero_to_one_based))    # id
        out.write('{} '.format(connection[1][0])) # x coord
        out.write('{} '.format(connection[1][1])) # y coord
        out.write('{} '.format(connection[1][2])) # z coord
        out.write('{}'.format(connection[2]))     # area
        out.write('\n')
    out.close()

if __name__ == "__main__":

    h5 = False
    if h5:
        sys.exit('HDF5 formatted explicit unstructured grids currently not '+
                 'supported by struct_grid_to_ugrid_explicit()')
        out_name = "mesh_uge.h5"
    else:
        out_name = "mesh.uge"

    n, d, origin = get_dimensions()
    struct_grid_to_ugrid_explicit(out_name, h5, n, d, origin)

    print('Done with unstructured explicit grid.')
