"""
Created on Thu Sep 10 12:04:03 2020

@author: rosie
"""

import os, sys
from pydfnworks import * 

src_path = os.getcwd()

jobname = f"{src_path}/output"   
        
#dfnFlow_file = f"{src_path}/dfn_explicit.in"
DFN = DFNWORKS(jobname)
#               dfnFlow_file=dfnFlow_file)

DFN.params['domainSize']['value'] = [1000.0, 1000.0, 1000.0]
DFN.params['h']['value'] = 5.0

DFN.add_user_fract(shape='ell',
                    radii=600,
                    translation=[-400, 0, 200],
                    normal_vector=[30, 15, 60],
                    number_of_vertices=5,
                    aperture=1.0e-3)

DFN.add_user_fract(shape='ell',
                    radii=1000,
                    translation=[0, 0, 0],
                    normal_vector=[95, 5, 0],
                    number_of_vertices=5,
                    aperture=1.0e-3)

DFN.add_user_fract(shape='ell',
                    radii=600,
                    aspect_ratio=1,
                    translation=[400, 0, 200],
                    normal_vector=[30, 15, 60],
                    number_of_vertices=5,
                    aperture=1.0e-3)

DFN.add_user_fract(shape='ell',
                    radii=600,
                    aspect_ratio=1,
                    translation=[400, 0, -400],
                    normal_vector=[30, 15, 60],
                    number_of_vertices=5,
                    aperture=5.0e-5)

DFN.make_working_directory(delete=True)
DFN.print_domain_parameters()
DFN.check_input()
DFN.create_network()
DFN.dump_hydraulic_values()
#DFN.mesh_network(uniform_mesh=True)

#DFN.dfn_flow()
