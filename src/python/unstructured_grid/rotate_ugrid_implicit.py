import sys
from h5py import *
import numpy
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
        rotation_radians = angleXZ/180*math.pi
        new_x = math.cos(rotation_radians)*x-math.sin(rotation_radians)*z
        new_z = math.sin(rotation_radians)*x+math.cos(rotation_radians)*z
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

filename = ''

if len(sys.argv) > 1:
  filename = sys.argv[1]

origin = [0.,0.,0.]
angleXY = 0.
angleXZ = 0.
angleYZ = 0.

h5file = File(filename,mode='r+')

vertices = h5file['Domain/Vertices']

for i in range(vertices.shape[0]):
    print(vertices[i])
    vertex = rotate_around_point(origin,angleXY,angleXZ,angleYZ,vertices[i])
    vertices[i] = vertex
    print(vertex)

h5file.close()
