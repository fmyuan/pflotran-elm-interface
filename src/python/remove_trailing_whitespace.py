import sys
import shutil
import os
import re

path = ''

def get_filename(root,suffix):
    filename = []
    filename.append(root)
    filename.append(suffix)
    filename = '.'.join(filename)
    return filename
  
def copy_files(source_file_roots):
    for root in source_file_roots:
      filename = get_filename(root,'F90')
      shutil.copy(filename,path+filename)
  
def open_f(filename):
    return open(path+filename,'r')
def open_f2(filename):
    return open(path+filename+'.tmp','w')
def mv_file(filename):
    # using shutil.move adds ^M to end of lines.
    shutil.copy(path+filename+'.tmp',path+filename)
    os.remove(path+filename+'.tmp')

def refactor_file(filename):
    white_space_removed = 0
    lines_shortened = 0
    f = open_f(filename)
    f2 = open_f2(filename)
    for line in f:
        len0 = len(line)
        line = line.rstrip()+'\n'
        len1 = len(line)
        if len0 > len1:
            lines_shortened += 1
            white_space_removed += len0-len1
        f2.write(line)
    f.close()
    f2.close()
    mv_file(filename)
    return white_space_removed, lines_shortened

# Obtain list of source files
source_file_roots = []
for line in open('pflotran_object_files.txt','r'):
    # find .o file
    # could use re.split() here, but too complicated.
    w = line.split('}')
    if len(w) == 2:
        w2 = w[1].split('.o')
#        print(w2[0])
        source_file_roots.append(w2[0])
source_file_roots.append('pflotran')
source_file_roots.append('pflotran_rxn')
source_file_roots.append('pflotran_derivative')

#clm-pflotran files
source_file_roots.append('../clm-pflotran/pflotran_model')
source_file_roots.append('../clm-pflotran/pflotran_interface_main')
source_file_roots.append('../clm-pflotran/mapping')
source_file_roots.append('../clm-pflotran/clm_pflotran_interface_data')
source_file_roots.append('unittests/test_characteristic_curves.pf')

#unit test files
source_file_roots.append('unittests/test_characteristic_curves_thermal.pf')
source_file_roots.append('unittests/test_eos_gas.pf')
source_file_roots.append('unittests/test_eos_water.pf')
source_file_roots.append('unittests/test_geometry.pf')
source_file_roots.append('unittests/test_material.pf')
source_file_roots.append('unittests/test_saturation_function.pf')
source_file_roots.append('unittests/test_string.pf')
source_file_roots.append('unittests/test_utility.pf')

# Alphabetize
source_file_roots.sort()
print(source_file_roots)
#os.remove('pflotran_provenance.F90')
file_count = 0
white_space_removed = 0
lines_shortened = 0
for root in source_file_roots:
    file_count += 1
    if not root.endswith('pf'):
        filename = get_filename(root,'F90')
    else:
        filename = root
    print(filename)
    n0, n1 = refactor_file(filename)
    white_space_removed += n0
    lines_shortened += n1

print('')
print('{} lines shortened.'.format(lines_shortened))
print('{} whitespace characters removed.'.format(white_space_removed))
print('done!')
