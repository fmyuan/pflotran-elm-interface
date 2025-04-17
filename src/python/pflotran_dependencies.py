# Author: Glenn Hammond
# Date: 06/03/13
# This python script calculates the dependencies between PFLOTRAN source files
# and writes the dependencies to pflotran_dependencies.txt.  The developer
# must then insert the contents of this file within the PFLOTRAN makefile.

# YOU MUST RUN THIS SCRIPT FROM WITHIN PFLOTRAN_DIR/src/pflotran
# > python ../python/pflotran_dependencies.py

import os
import sys
import shutil

run_test = False
if len(sys.argv) > 1:
    run_test = 'test' in sys.argv[1:]

try:
  pflotran_dir = os.environ['PFLOTRAN_DIR']
except KeyError:
  print('ERROR: PFLOTRAN_DIR must point to PFLOTRAN installation directory and be defined in system environment variables.')
  sys.exit(1)
sys.path.append(pflotran_dir + '/src/python/util')

from source_files import get_source_file_roots

def get_filename(root,suffix):
  filename = []
  filename.append(root)
  filename.append(suffix)
  filename = '.'.join(filename)
  return filename

# generate pflotran_provenance.F90 if it does not exit
filename = 'pflotran_provenance.F90'
if not os.path.exists(filename):
  shutil.copy('pflotran_no_provenance.F90',filename)

pflotran_rxn_list = []
remove_file_list = []
remove_file_list.append('logging')
pflotran_rxn_list.append(('constraint',remove_file_list))
remove_file_list = []
differing_pflotran_rxn_dependencies = dict(pflotran_rxn_list)

module_skip_list = ('hdf5','h5lt','petsc','clm_pflotran_interface_data', \
                    'ieee_arithmetic','iso_c_binding')

# Obtain list of source files
source_file_roots = get_source_file_roots()

# Alphabetize
source_file_roots.sort()
#print(source_file_roots)
f = open('pflotran_source_files.txt','w')
f.write('Filename                                          source   blank comment\n')
f.write('------------------------------------------------------------------------\n')
file_count = 0 
line_count = 0
blank_line_count = 0
comment_line_count = 0
for root in source_file_roots:
  file_count += 1
  f.write(get_filename(root,'F90'))
  line_count_in_file = 0
  blank_line_count_in_file = 0
  comment_line_count_in_file = 0
  filename = get_filename(root,'F90')
  try:
    for line in open(filename):
      stripped = line.lstrip()
      if stripped.startswith('!'):
        comment_line_count += 1
        comment_line_count_in_file += 1
      elif len(stripped) < 1:
        blank_line_count += 1
        blank_line_count_in_file += 1
      else:
        line_count += 1
        line_count_in_file += 1
  except Exception as e:
    print(f'ERROR: Unexpected error {e}')
    sys.exit(1)
  for i in range(50-len(filename)):
    f.write(' ')
  f.write('%6d  %6d  %6d\n' % (line_count_in_file,blank_line_count_in_file,
                            comment_line_count_in_file))
f.write('------------------------------------------------------------------------\n')
f.write('Number of files: %d\n' % file_count)
f.write('Total Line count: %d\n' % 
        (line_count+comment_line_count+blank_line_count))
f.write('Fortran source line count: %d\n' % line_count)
f.write('Fortran comment line count: %d\n' % comment_line_count)
f.write('Blank line count: %d\n' % blank_line_count)
f.write('Average source lines per file: %.1f\n' % 
        (float(line_count) / file_count))
f.close()

# Obtain a list of modules
module_list = []
for root in source_file_roots:
#  print(root)
  for line in open(get_filename(root,'F90')):
    if line.lstrip().startswith('module'):
      w = line.split()
      # skip 'module procedure'
      if not w[1].startswith('procedure'):
#        print('  '+w[1])
        module_link = []
        module_link.append(w[1])
        module_link.append(get_filename(root,'o'))
        module_list.append(module_link)
module_dictionary = dict(module_list)
#print(module_dictionary.keys())

dependency_filename_prefix = 'pflotran_dependencies'
f = open(dependency_filename_prefix+'.tmp','w')
# now loop over all source files and create the dependency list
for root in source_file_roots:
#  print(root)
  try:
    modules_to_remove = differing_pflotran_rxn_dependencies[root]
    num_times_to_print = 2
    f.write('ifdef PFLOTRAN_RXN_FLAG\n')
  except:
    num_times_to_print = 1
  for iprint in range(num_times_to_print):
    module_list = []
    for line in open(get_filename(root,'F90')):
      if line.lstrip().startswith('use '):
        w = line.split()
        # split again as there may be a ', only ...' after the module name
        # cannot use string.strip(',') as it will not work in the case of
        # XXX_module,only...
        w = w[1].split(',')
        # skip modules
        if not w[0].startswith(module_skip_list):
          module_list.append(w[0])
    # remove duplicate modules
    module_list = set(module_list)
    file_list = []
    for module in module_list:
      try:
        key = module_dictionary[module]
      except:
        # need to skip hdf5
        if not module.startswith(module_skip_list):
          print('ERROR: Module "%s" not found in dictionary.\n' % module)
          print(root, module)
          sys.exit(1)
#      print(key)
      file_list.append(key)
    # remove duplicates first as it will destroy an sorting
    file_set = set(file_list)
    # convert back to list
    file_list = list(file_set)
    # sort
    sorted_file_list = sorted(file_list)
    # remove the root file if it is listed
    filename = get_filename(root,'o')
    if filename in sorted_file_list:
      sorted_file_list.remove(filename)
      print('Removing %s from own dependency.\n' % filename)
#    print(sorted_file_list)
    # remove files that are not needed for pflotran rxn
    if num_times_to_print == 2 and iprint == 0:
      for root2 in modules_to_remove:
        sorted_file_list.remove(get_filename(root2,'o'))
#      print(sorted_file_list)
    multiple_files_per_line = False
    if multiple_files_per_line:
      string = '%s :' % get_filename(root,'o')
      f.write('%s' % string)
      string_len = len(string)
      indentation = string_len
      for file in sorted_file_list:
        string = ' %s' % file
        word_len = len(string)
        string_len += word_len
        if string_len > 78:
          f.write(' \\\n')
          for i in range(indentation):
            f.write(' ')
          string_len  = indentation + word_len
        f.write('%s' % string)
      f.write('\n')
    else:
      num_files = len(sorted_file_list)
      if num_files > 0:
        f.write('%s : \\\n' % get_filename(root,'o'))
      else:
        f.write('%s :\n' % get_filename(root,'o'))
      for ifile in range(num_files):
        f.write('  %s' % sorted_file_list[ifile])
        if ifile < num_files-1:
          f.write(' \\')
        f.write('\n')
    if num_times_to_print == 2 and iprint == 0:
      f.write('else\n')
    elif num_times_to_print == 2: 
      f.write('endif\n')
f.close()  

if not run_test:
  shutil.move(dependency_filename_prefix+'.tmp',
              dependency_filename_prefix+'.txt')

print('done!')
