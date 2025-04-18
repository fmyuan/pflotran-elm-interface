import sys
import shutil
import os
import re

try:
  pflotran_dir = os.environ['PFLOTRAN_DIR']
except KeyError:
  print('PFLOTRAN_DIR must point to PFLOTRAN installation directory and be defined in system environment variables.')
  sys.exit(1)
sys.path.append(pflotran_dir + '/src/python/util')

from source_files import *

path = ''

def get_filename(root,suffix):
    filename = []
    filename.append(root)
    filename.append(suffix)
    filename = '.'.join(filename)
    return filename
  
def open_f(filename):
    return open(path+filename,'r')

def check_file_line_length(filename,max_line_length):
    lines_greater_than_max_length = 0
    f = open_f(filename)
    iline = 0
    for line in f:
        iline += 1
        if len(line.split('!')[0]) > max_line_length:
            lines_greater_than_max_length += 1
            print_error(filename,line,iline)
    f.close()
    return lines_greater_than_max_length

def print_error(filename,line,iline):
        print(f'{filename}({iline}) : {len(line)} characters')

def main():

    max_line_length = 132
    source_file_roots = get_source_file_roots() + \
                        get_clm_pflotran_source_file_roots() + \
                        get_unit_test_files()

    # remove pflotran_provenance.o is it does not exit before compilation
    source_file_roots.remove('pflotran_provenance')

    total_num_lines_failing = 0
    num_files_impacted = 0
    violating_filenames = []
    violating_file_line_counts = []
    for root in source_file_roots:
        if not root.endswith('pf'):
            filename = get_filename(root,'F90')
        else:
            filename = root
        num_lines_failing = check_file_line_length(filename,max_line_length)
        total_num_lines_failing += num_lines_failing
        print(f'{filename}: {num_lines_failing}')
        if num_lines_failing > 0:
            violating_filenames.append(filename)
            violating_file_line_counts.append(num_lines_failing)
            num_files_impacted += 1

    if num_files_impacted > 0:
        print(f'\n{num_files_impacted} file(s) violating the maximum ' + \
              'compilable (excluding comments) line length of\n' +
              f'{max_line_length} characters:\n')
        for i in range(num_files_impacted):
            print(f'    {violating_filenames[i]} : {violating_file_line_counts[i]}')
        print()
    else: 
        print('\nNo long lines found.\n')

if __name__ == "__main__":
    try:
        suite_status = main()
        print("success running "+__file__)
        sys.exit(suite_status)
    except Exception as error:
        print(str(error))
#        if cmdl_options.backtrace:
#            traceback.print_exc()
#        traceback.print_exc()
        print("failure running "+__file__)
        sys.exit(1)

    print('done!')
