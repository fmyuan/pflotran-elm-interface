import os
import sys
try:
    pflotran_dir = os.environ['PFLOTRAN_DIR']
except KeyError:
    print('PFLOTRAN_DIR must point to PFLOTRAN installation directory and be defined in system environment variables.')
    sys.exit(1)
sys.path.append(pflotran_dir + '/src/python/util')
from source_files import *

import shutil
import fnmatch
import re

myflags = re.I

def mv_file(old_filename,new_filename):
    # using shutil.move adds ^M to end of lines.
    shutil.copy(old_filename,new_filename)
    os.remove(old_filename)

regex_string = r'(?P<pre>\s*)(?P<variable>\S*)\s*=\s*' + \
               r'(?P<null_expression>PETSC_NULL_\w*)(?P<ending>.*)'
regex_ = re.compile(regex_string,flags=myflags)

def refactor(filename,overwrite):
    fin = open(filename,'r')
    fout = open(filename+'.tmp','w')
    prev_line = ''
    iline = 0
    for line in fin:
        iline += 1
        m = re.search(regex_,line)
        if m:
            if m.group('null_expression').find('CHARACTER') == -1:
                line = m.group('pre') + \
                   'PetscObjectNullify(' + \
                   m.group('variable') + \
                   ')' + \
                   m.group('ending') + \
                   '\n'
#            print(line)
        fout.write(line)
    fin.close()
    fout.close()
    if overwrite:
        mv_file(filename+'.tmp',filename)

def main():
    overwrite = False
    overwrite = True
    test = False
#    test = True
    if test:
        filename_list = ['test']
    else:
        filename_list = get_source_file_roots()
        filename_list.remove('petsc_utility')
    print(filename_list)
    for filename in filename_list:
        filename += '.F90'
        print(f'refactoring {filename}')
        refactor(filename,overwrite)

if __name__ == "__main__":
    try:
        suite_status = main()
        print("success")
        sys.exit(suite_status)
    except Exception as error:
        print(str(error))
#        if cmdl_options.backtrace:
#            traceback.print_exc()
#        traceback.print_exc()
        print("failure")
        sys.exit(1)

print('done')
