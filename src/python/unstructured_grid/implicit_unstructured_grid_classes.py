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
        
class Grid:
    
    def __init__(self,filename):
        self.filename = filename
        self.cells = []
        self.vertices = []
        
    def get_vertices(self):
        return self.vertices
        
    def get_cells(self):
        return self.cells
        
    def read_grid(self):
        print(self.filename)
        if self.filename.endswith('.h5'):
            h5file = h5py.File(self.filename,'r')
        else:
            print_err_msg('Only HDF5 files are currently supported')
            
        cells = h5file['Domain/Cells']
        cells = np.array(cells[:,:])        
        for i in range(cells.shape[0]):
            if i % 1000000 == 0:
                print('cells: ',i)
            self.cells.append(Cell(i+1,cells[i,0],cells[i,1:]))

        vertices = h5file['Domain/Vertices']
        vertices = np.array(vertices[:,:])        
        for i in range(vertices.shape[0]):
            if i % 1000000 == 0:
                print('vertices: ',i)
            self.vertices.append(Vertex(i+1,vertices[i,0],vertices[i,1],
                                        vertices[i,2]))
        h5file.close()
        print('Finished reading grid')
        
    def cross_reference(self):
        print('Cross Referencing Grid')
        i = 0
        for cell in self.cells:
            if i % 1000000 == 0:
                print('xref grid: ',i)
            for vertex_id in cell.vertex_ids:
                self.vertices[vertex_id-1].add_cell(cell.get_id())
            i += 1
                
    def print_cross_reference(self,fid):
        for vertex in self.vertices:
            vertex.print_cells(fid)
            
def read_regions_from_file(filename):
    
    if filename.endswith('.h5'):
        h5file = h5py.File(filename,'r')
    else:
        print_err_msg('Only HDF5 files are currently supported for regions')
            
    sideset_list = []
    regions_group = h5file['Regions']
    region_names = regions_group.keys()
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
        
    def cross_reference(self,cell_list,vertex_list):
        icount = 0
        for face in self.faces:
            icount += 1
            if icount % 100 == 0:
                print('xref face: ',icount)
            list_of_cell_lists = []
            for vertex_id in face.get_vertex_ids():
                vertex = vertex_list[vertex_id-1]
                list_of_cell_lists.append(vertex.get_cell_ids())
            in_all = set.intersection(*[set(x) for x in list_of_cell_lists])
            face.cell_ids = in_all
                
    def print_faces(self,fid,cell_list):
        write(fid,'Sideset: {}\n'.format(self.name))
        for face in self.faces:
            face.print_face(fid,cell_list)

class Cell:
    
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
            print_err_msg('Unknown cell type in cell {}: itype = {}'.
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
    
    def get_id(self):
        return self.id

    def get_vertex_ids(self):
        return self.vertex_ids

    def get_ctype(self):
        return self.ctype
    
    def print_cell(self,fid):
        write(fid,'  Cell {} [{}]:'.
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
        self.cell_ids = []
        
    def get_id(self):
        return self.id

    def get_vertex_ids(self):
        return self.vertex_ids

    def get_cell_ids(self):
        return self.cell_ids

    def get_ctype(self):
        return self.ctype    

    def print_face(self,fid,cell_list):  
        vertex_ids = self.get_vertex_ids()
        cell_ids = self.get_cell_ids()
        write(fid,'Face {} [{}] : {} vertices, {} cell(s)\n'.
              format(self.get_id(),self.get_ctype(),len(vertex_ids),
                     len(cell_ids)))
        write(fid,'  Vertices :')
        for vertex_id in vertex_ids:
            write(fid,' {}'.format(vertex_id))
        write(fid,'\n')
        for cell_id in cell_ids:
            cell_list[cell_id-1].print_cell(fid)

    
class Vertex:

    def __init__(self,id,x,y,z):
        self.id = id        
        self.x = x
        self.y = y
        self.z = z
        self.cell_ids = []
        
    def get_id(self):
        return self.id
        
    def add_cell(self,cell_id):
        self.cell_ids.append(cell_id)
        
    def get_cell_ids(self):
        return self.cell_ids
    
    def print_cells(self,fid):
        write(fid,'Vertex {} :'.format(self.get_id()))
        for cell_id in self.cell_ids:
            write(fid,' {}'.format(cell_id))
        write(fid,'\n')
    

