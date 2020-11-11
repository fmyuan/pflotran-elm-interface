"""
Description: Classes defining cells, faces, and vertices for an implicit 
             unstructured grid.
Author: Glenn Hammond
"""

import sys
import os

try:
    pflotran_dir = os.environ['PFLOTRAN_DIR']
except KeyError:
    print('PFLOTRAN_DIR must point to PFLOTRAN installation directory and '+
          'be defined in system environment variables.')
    sys.exit(1)
sys.path.append(pflotran_dir + '/src/python')

from common.progress_bar import ProgressBar
from common.output import Output
from common.error_messaging import print_err_msg

import h5py
import numpy as np
        
class Grid:
    
    def __init__(self,filename):
        self.filename = filename
        self.cells = []
        self.vertices = []
        
    def get_vertices(self):
        return self.vertices
        
    def get_cells(self):
        return self.cells
    
    def get_num_cells(self):
        return len(self.cells)
        
    def read_grid(self,output):
        
        output.print_and_return('Reading grid from {}'.format(self.filename))
        if self.filename.endswith('.h5'):
            h5file = h5py.File(self.filename,'r')
        else:
            print_err_msg('Only HDF5 files are currently supported')
            
        output.print_and_return('Reading cells')
        cells = h5file['Domain/Cells']
        cells = np.array(cells[:,:])   
        pb = ProgressBar(cells.shape[0])
        for i in range(cells.shape[0]):
            pb.increment()
            self.cells.append(Cell(i+1,cells[i,0],cells[i,1:]))

        output.print_and_return('Reading vertices')
        vertices = h5file['Domain/Vertices']
        vertices = np.array(vertices[:,:])   
        pb = ProgressBar(vertices.shape[0])
        for i in range(vertices.shape[0]):
            pb.increment()
            self.vertices.append(Vertex(i+1,vertices[i,0],vertices[i,1],
                                        vertices[i,2]))
        h5file.close()
        output.print_and_return('Finished reading grid')
        
    def cross_reference(self,output):
        output.print_and_return('Cross Referencing Grid')
        pb = ProgressBar(len(self.cells))
        for cell in self.cells:
            pb.increment()
            for vertex_id in cell.vertex_ids:
                self.vertices[vertex_id-1].add_cell(cell.get_id())
                
    def print_cross_reference(self,output):
        output.print_and_return('Printing vertex-cell cross reference')
        for vertex in self.vertices:
            vertex.print_cells(output)
            
def read_regions_from_file(output,filename):
    
    if filename.endswith('.h5'):
        h5file = h5py.File(filename,'r')
    else:
        print_err_msg('Only HDF5 files are currently supported for regions')
            
    output.print_and_return('Reading regions with vertex ids from {}'.format(filename))
    sideset_list = []
    try:
        regions_group = h5file['Regions']
    except:
        output.print_and_return('{} lacks a Regions group'.format(filename))
        return []
    region_names = regions_group.keys()
    for name in region_names:
        try:
            temp_name = 'Regions/'+name+'/Vertex Ids'
            vertices = h5file[temp_name]
        except:
            output.print_and_return('Region "{}" lacks vertex ids'.format(name))
            continue
        output.print_and_return('Reading vertex ids from region: {}'.format(name))
        sideset = SideSet(name)
        sideset.parse_sideset(vertices)
        sideset_list.append(sideset)

    h5file.close()
    return sideset_list
          
class SideSet:
    
    def __init__(self,name):
        self.name = name
        self.faces = []
        self.cell_ids = []
        
    def parse_sideset(self,vertices):
        vertices = np.array(vertices[:,:])        
        for i in range(vertices.shape[0]):
            self.faces.append(Face(i+1,vertices[i,0],vertices[i,1:]))
        
    def cross_reference(self,output,cell_list,vertex_list):
        output.print_and_return('Cross referencing sideset: {}'.format(self.name))
        pb = ProgressBar(len(self.faces))
        for face in self.faces:
            pb.increment()
            list_of_cell_lists = []
            for vertex_id in face.get_vertex_ids():
                vertex = vertex_list[vertex_id-1]
                list_of_cell_lists.append(vertex.get_cell_ids())
            in_all = set.intersection(*[set(x) for x in list_of_cell_lists])
            face.cell_ids = in_all
        output.print_and_return()

    def print_info(self,output):
        output.print_and_return('Sideset: {}'.format(self.name))
        output.print_and_return('  Face -> Cell Counts: ')
        counts = self.get_face_cell_counts()
        output.print_and_return('     total # faces = {}'.format(counts[0]))
        output.print_and_return('     with  0 cells = {}'.format(counts[1]))
        output.print_and_return('     with  1 cell  = {}'.format(counts[2]))
        output.print_and_return('     with >1 cell  = {}'.format(counts[3]))
                
    def print_faces(self,output,cell_list):
        output.set_indentation(2)
        output.print_and_return('Faces:'.format(self.name))
        output.set_indentation(4)
        for face in self.faces:
            face.print_face(output,cell_list)
        output.clear_indentation()

    def get_face_cell_counts(self):
        icount = [0]*4
        icount[0] = len(self.faces)
        for face in self.faces:
            num_cells = len(face.cell_ids)
            if num_cells == 0:
                icount[1] += 1
            elif num_cells == 1:
                icount[2] += 1
            else:
                icount[3] += 1
        return icount

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
    
    def print_cell(self,output):
        strings = ['Cell {} [{}]:'.format(self.get_ctype(),self.get_id())]
        for vertex_id in self.get_vertex_ids():
            strings.append(' {}'.format(vertex_id))
        output.print_and_return(''.join(strings))
    
    
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

    def print_face(self,output,cell_list):  
        vertex_ids = self.get_vertex_ids()
        cell_ids = self.get_cell_ids()
        output.print_and_return('Face {} [{}] : {} vertices, {} cell(s)'.
              format(self.get_id(),self.get_ctype(),len(vertex_ids),
                     len(cell_ids)))
        output.add_indentation(2)
        strings = ['Vertices :']
        for vertex_id in vertex_ids:
            strings.append(' {}'.format(vertex_id))
        output.print_and_return(''.join(strings))
        for cell_id in cell_ids:
            cell_list[cell_id-1].print_cell(output)
        output.subtract_indentation(2)

    
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
    
    def print_cells(self,output):
        strings = ['Vertex {} :'.format(self.get_id())]
        for cell_id in self.cell_ids:
            strings.append(' {}'.format(cell_id))
        output.print_and_return(''.join(strings))
    

