"""
Description: Reads an implicit unstructured grid and ensures that associated
             regions are properly set up.
Author: Glenn Hammond
"""
import sys
import os
import traceback

try:
  pflotran_dir = os.environ['PFLOTRAN_DIR']
except KeyError:
  print('PFLOTRAN_DIR must point to PFLOTRAN installation directory and be defined in system environment variables.')
  sys.exit(1)
sys.path.append(pflotran_dir + '/src/python')

from unstructured_grid.implicit_unstructured_grid_classes import *

test = True

def main():
    
    fid = None
    filename = 'parse.stdout'
    fid = open(filename,'w')
    
    status = 0
    if test:
        filename = 'mixed.h5'
    else:
        filename = 'Thorne_mesh_10m-matID-river.h5'
    grid = Mesh(filename)
    grid.read_mesh()
    grid.cross_reference()
    if test:
        grid.print_cross_reference(fid)
    
    if test:
        filename = 'regions.h5'
    sideset_list = read_regions_from_file(filename,test)
    for sideset in sideset_list:
        sideset.cross_reference(grid.get_elements(),grid.get_vertices())
        sideset.print_faces(fid,grid.get_elements())
    
    if fid:
        fid.close()
    return status
    
if __name__ == "__main__":
    try:
        status = main()
        print("success")
    except Exception as error:
        print(str(error))
#        if cmdl_options.backtrace:
#            traceback.print_exc()
        traceback.print_exc()
        print("failure")
        sys.exit(1)
