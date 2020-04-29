import sys
import h5py
import numpy as np
import math

def rotate(angleXY,angleXZ,angleYZ,vertex):
    num_angle = 0
    if abs(angleXY) > 0.:
        num_angle += 1
    if abs(angleXZ) > 0.:
        num_angle += 1
    if abs(angleYZ) > 0.:
        num_angle += 1
    if num_angle > 1:
        sys.exit('Too many angles')
    x,y,z = vertex
    if abs(angleXY) > 0.:
        rotation_radians = angleXY/180*math.pi
        new_x = math.cos(rotation_radians)*x-math.sin(rotation_radians)*y
        new_y = math.sin(rotation_radians)*x+math.cos(rotation_radians)*y
        x = new_x
        y = new_y
    if abs(angleXZ) > 0.:
        #Moise: changed to make rotation according to right hand rule around Oy
        rotation_radians = angleXZ/180*math.pi
        new_x = math.cos(rotation_radians)*x+math.sin(rotation_radians)*z
        new_z = -math.sin(rotation_radians)*x+math.cos(rotation_radians)*z
        x = new_x
        z = new_z
    if abs(angleYZ) > 0.:
        rotation_radians = angleYZ/180*math.pi
        new_y = math.cos(rotation_radians)*y-math.sin(rotation_radians)*z
        new_z = math.sin(rotation_radians)*y+math.cos(rotation_radians)*z
        y = new_y
        z = new_z
    return [x,y,z]
    
def translate(dx,dy,dz,vertex):
    x,y,z = vertex
    vertex[0] += dx
    vertex[1] += dy
    vertex[2] += dz
    return [x,y,z]
    
def rotate_around_point(point,angleXY,angleXZ,angleYZ,vertex):
    vertex = translate(-point[0],-point[1],-point[2],vertex)
    vertex = rotate(angleXY,angleXZ,angleYZ,vertex)
    vertex = translate(point[0],point[1],point[2],vertex)
    return vertex
    

  
def rotate_ugrid(input_grid, output_grid, angle, origin):
    """
    Rotate the input grid and store it into the output file
    Rotate one angle per angle only, since the rotation order matters
    If the output file is the same as the input, the input grid is erased
    input_grid: input grid file
    output_grid: output grid file
    angle = (angleXY,angleXZ,angleYZ)
    origin = (x,y,z)
    """
    
    if input_grid != output_grid:
        src = h5py.File(input_grid, mode='r')
        out = h5py.File(output_grid, mode='w')
        out.create_dataset('Domain/Cells', data=src['Domain/Cells'])
        out.create_dataset('Domain/Vertices', data=src['Domain/Vertices'])
        src.close()
        out.close()
    
    src = h5py.File(output_grid,mode='r+')
    
    angleXY, angleXZ, angleYZ = angle
    
    vertices = np.copy(src['Domain/Vertices'])
    for i in range(vertices.shape[0]):
        vertex = rotate_around_point(origin,angleXY,angleXZ,angleYZ,vertices[i])
        vertices[i] = vertex
        
    del src['Domain/Vertices']
    src.create_dataset('Domain/Vertices', data=vertices)

    src.close()
    
    return

    
def compute_perm_tensor_component(K,angle):
    """
    Compute the tensor component to express the kx, ky, kz permeabilty of the material in the simualtion base
    """
    angleXY,angleXZ,angleYZ = angle
    num_angle = 0
    if abs(angleXY) > 0.:
        num_angle += 1
    if abs(angleXZ) > 0.:
        num_angle += 1
    if abs(angleYZ) > 0.:
        num_angle += 1
    if num_angle > 1:
        sys.exit('Too many angles')
    
    kxx,kxy,kxz,kyy,kyz,kzz = K
    
    Kmat = np.array([[kxx, kxy, kxz], [kxy, kyy, kyz], [kxz,kyz, kzz]])
    
    if abs(angleXY) > 0.:
        angleXY = angleXY/180*np.pi
        P = [[np.cos(angleXY),-np.sin(angleXY),0], [np.sin(angleXY), np.cos(angleXY), 0], [0,0,1]]
    elif abs(angleXZ) > 0.:
        angleXZ = angleXZ/180*np.pi
        P = [[np.cos(angleXZ), 0, np.sin(angleXZ)], [0, 1, 0], [-np.sin(angleXZ), 0, np.cos(angleXZ)]]
    elif abs(angleYZ) > 0.:
        angleYZ = angleYZ/180*np.pi
        P = [[1, 0, 0], [0, np.cos(angleYZ), -np.sin(angleYZ)], [0, np.sin(angleYZ), np.cos(angleYZ)]]
    else:
        P = [[1, 0, 0], [0,1,0], [0,0,1]]
        
    res = np.dot(P,Kmat)
    res = np.dot(res,np.transpose(P))
    
    K = (res[0,0], res[0,1], res[0,2], res[1,1], res[1,2], res[2,2])
    
    return K


if __name__ == "__main__":
    input_grid = "mesh_ref_ugi.h5"
    output_grid = "mesh_rotated_ugi.h5"
    angleXY = 30.
    angleXZ = 0.
    angleYZ = 0.
    angle = (angleXY, angleXZ, angleYZ)
    origin = (0,0,0)
    rotate_ugrid(input_grid, output_grid, angle, origin)
    K = (1e-5,0,0,1e-6,0,3e-6)
    compute_perm_tensor_component(K,angle)
    print(K)





