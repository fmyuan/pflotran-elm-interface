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

def refactor_file(filename,overwrite_file):
    tabs_removed = 0
    lines_converted = 0
    f = open_f(filename)
    if overwrite_file:
        f2 = open_f2(filename)
    for line in f:
        line2 = line.replace('\t','  ')
        if not line2 == line:
            lines_converted += 1
            for character in line:
                if character == '\t':
                    tabs_removed += 1
        if overwrite_file:
            f2.write(line2)
    f.close()
    if overwrite_file:
        f2.close()
        mv_file(filename)
    return tabs_removed, lines_converted

    # Alphabetize
    source_file_roots.sort()
    #print(source_file_roots)
    return source_file_roots

def main():

    overwrite_file = False
    if len(sys.argv) > 1:
        args = []
        num_recognized = 0
        for arg in sys.argv[1:]:
            args.append(arg.strip())
        if 'overwrite' in args:
            num_recognized += 1
            overwrite_file = True
        if not num_recognized == len(args):
            print(args)
            raise Exception('Unrecognized arguments.')
     
    source_file_roots = get_source_file_roots() + \
                        get_clm_pflotran_source_file_roots() + \
                        get_unit_test_files()
    
    # remove pflotran_provenance.o is it does not exit before compilation
    source_file_roots.remove('pflotran_provenance')
        
    file_count = 0
    tabs_removed = 0
    lines_converted = 0
    num_files_impacted = 0
    refactored_files = []
    for root in source_file_roots:
        file_count += 1
        if not root.endswith('pf'):
            filename = get_filename(root,'F90')
        else:
            filename = root
        n0, n1 = refactor_file(filename,overwrite_file)
        if n0 > 0:
            refactored_files.append(filename)
        print('{}: {} {}'.format(filename,n0,n1))
        if n0 > 0:
            num_files_impacted += 1
        tabs_removed += n0
        lines_converted += n1

    if tabs_removed > 0:
        if len(refactored_files) > 1:
            string = ', '.join(refactored_files)
            string = 'files impacted: ' + string
        else:
            string = refactored_files[0]
            string = 'file impacted: ' + string
        print('\n{} {}'.format(num_files_impacted,string))
        if overwrite_file:
            print('{} lines refactored'.format(lines_converted))
            print('{} tab characters removed\n'.format(tabs_removed))
        else:
            print('{} lines impacted'.format(lines_converted))
            print('{} tab characters found\n'.format(tabs_removed))
    else:
        print('\nNo tab characters found.\n')

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
