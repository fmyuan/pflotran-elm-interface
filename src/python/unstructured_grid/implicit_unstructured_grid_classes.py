"""
Description: Classes defining cells, faces, and vertices for an implicit 
             unstructured grid.
Author: Glenn Hammond
"""

import sys
import h5py
import numpy as np

def print_err_msg(*strings):
    list = []
    for string in strings:
        list.append(string)
    string = ''.join(list)
    if True:
        raise Exception(string)
    else:
        print(string)
        sys.exit(1)
        
def write(fid,string):
    if fid:
        fid.write(string)
    sys.stdout.write(string)
        
class Mesh:
    
    def __init__(self,filename):
        self.filename = filename
        self.elements = []
        self.vertices = []
        
    def get_vertices(self):
        return self.vertices
        
    def get_elements(self):
        return self.elements
        
    def read_mesh(self):
        print(self.filename)
        if self.filename.endswith('.h5'):
            h5file = h5py.File(self.filename,'r')
        else:
            print_err_msg('Only HDF5 files are currently supported')
            
        cells = h5file['Domain/Cells']
        cells = np.array(cells[:,:])        
        for i in range(cells.shape[0]):
            if i % 1000000 == 0:
                print('elements: ',i)
            self.elements.append(Element(i+1,cells[i,0],cells[i,1:]))

        vertices = h5file['Domain/Vertices']
        vertices = np.array(vertices[:,:])        
        for i in range(vertices.shape[0]):
            if i % 1000000 == 0:
                print('vertices: ',i)
            self.vertices.append(Vertex(i+1,vertices[i,0],vertices[i,1],
                                        vertices[i,2]))
        h5file.close()
        print('Mesh read')
        
    def cross_reference(self):
        print('Cross Referencing Mesh')
        i = 0
        for element in self.elements:
            if i % 1000000 == 0:
                print('xref mesh: ',i)
            for vertex_id in element.vertex_ids:
                self.vertices[vertex_id-1].add_element(element.get_id())
            i += 1
                
    def print_cross_reference(self,fid):
        for vertex in self.vertices:
            vertex.print_elements(fid)
            
def read_regions_from_file(filename,test):
    
    if filename.endswith('.h5'):
        h5file = h5py.File(filename,'r')
    else:
        print_err_msg('Only HDF5 files are currently supported for regions')
            
    sideset_list = []
    regions_group = h5file['Regions']
    region_names = regions_group.keys()
    if not test:
        region_names = ['Mass1_35']
    for name in region_names:
        try:
            temp_name = 'Regions/'+name+'/Vertex Ids'
            print(temp_name)
            vertices = h5file[temp_name]
        except:
            print('Region "{}" lacks vertex ids'.format(name))
            continue
        print('Creating sideset: ',name)
        sideset = SideSet(name)
        sideset.parse_sideset(vertices)
        sideset_list.append(sideset)

    h5file.close()
    return sideset_list
          
class SideSet:
    
    def __init__(self,name):
        self.name = name
        self.faces = []
        
    def parse_sideset(self,vertices):
        vertices = np.array(vertices[:,:])        
        for i in range(vertices.shape[0]):
            self.faces.append(Face(i+1,vertices[i,0],vertices[i,1:]))
        print(self.faces[0].get_vertex_ids())
        
    def cross_reference(self,element_list,vertex_list):
        icount = 0
        for face in self.faces:
            icount += 1
            if icount % 100 == 0:
                print('xref face: ',icount)
            list_of_element_lists = []
            for vertex_id in face.get_vertex_ids():
                vertex = vertex_list[vertex_id-1]
                list_of_element_lists.append(vertex.get_element_ids())
            in_all = set.intersection(*[set(x) for x in list_of_element_lists])
            face.element_ids = in_all
                
    def print_faces(self,fid,element_list):
        write(fid,'Sideset: {}\n'.format(self.name))
        for face in self.faces:
            face.print_face(fid,element_list)

class Element:
    
    def __init__(self,id,itype,vertex_ids):
        self.id = id
        self.itype = itype
        if itype == 4:
            self.ctype = 'T'
        elif itype == 5:
            self.ctype = 'P'
        elif itype == 6:
            self.ctype = 'W'
        elif itype == 8:
            self.ctype = 'H'
        else:
            print_err_msg('Unknown element type in element {}: itype = {}'.
                          format(self.id,self.itype))
        # only store non-zero vertex ids
        i = 0
        for vertex_id in vertex_ids:
            if vertex_id == 0:
                break
            i += 1
        self.vertex_ids = np.array(vertex_ids[:i])
        
    def get_ctype(self):
        return self.ctype
    
    def get_index(self):
        return self.id-1
    
    def get_id(self):
        return self.id

    def get_vertex_ids(self):
        return self.vertex_ids

    def get_ctype(self):
        return self.ctype
    
    def print_element(self,fid):
        write(fid,'  Element {} [{}]:'.
                         format(self.get_ctype(),self.get_id()))
        for vertex_id in self.get_vertex_ids():
            write(fid,' {}'.format(vertex_id))
        write(fid,'\n')
        
    
class Face:
    
    def __init__(self,id,itype,vertex_ids):
        self.id = id
        self.itype = itype
        if itype == 3:
            self.ctype = 'T'
        elif itype == 4:
            self.ctype = 'Q'
        else:
            print_err_msg('Unknown face type in face {}: itype = {}'.
                          format(self.id,self.itype))
        # only store non-zero vertex ids
        i = 0
        for vertex_id in vertex_ids:
            if vertex_id == 0:
                break
            i += 1
        self.vertex_ids = np.array(vertex_ids[:i])
        self.element_ids = []
        
    def get_index(self):
        return self.id-1
    
    def get_id(self):
        return self.id

    def get_vertex_ids(self):
        return self.vertex_ids

    def get_element_ids(self):
        return self.element_ids

    def get_ctype(self):
        return self.ctype    

    def print_face(self,fid,element_list):  
        vertex_ids = self.get_vertex_ids()
        element_ids = self.get_element_ids()
        write(fid,'Face {} [{}] : {} vertices, {} element(s)\n'.
              format(self.get_id(),self.get_ctype(),len(vertex_ids),
                     len(element_ids)))
        write(fid,'  Vertices :')
        for vertex_id in vertex_ids:
            write(fid,' {}'.format(vertex_id))
        write(fid,'\n')
        for element_id in element_ids:
            element_list[element_id-1].print_element(fid)

    
class Vertex:

    def __init__(self,id,x,y,z):
        self.id = id        
        self.x = x
        self.y = y
        self.z = z
        self.element_ids = []
        
    def get_index(self):
        return self.id-1
        
    def get_id(self):
        return self.id
        
    def add_element(self,element_id):
        self.element_ids.append(element_id)
        
    def get_element_ids(self):
        return self.element_ids
    
    def print_elements(self,fid):
        write(fid,'Vertex {} :'.format(self.get_id()))
        for element_id in self.element_ids:
            write(fid,' {}'.format(element_id))
        write(fid,'\n')
    

