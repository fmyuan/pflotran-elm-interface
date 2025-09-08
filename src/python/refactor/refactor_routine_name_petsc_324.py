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

routine_list = []
routine_list.append(['MatSetValuesLocal','PUMSetValuesLocal'])
routine_list.append(['MatSetValuesBlockedLocal','PUMSetValuesBlockedLocal'])
routine_list.append(['MatSetValues','PUMSetValues'])
routine_list.append(['MatSetValue','PUMSetValue'])
routine_list.append(['MatSetValuesBlocked','PUMSetValuesBlocked'])
routine_list.append(['MatGetValuesLocal','PUMGetValuesLocal'])
routine_list.append(['MatSeqAIJGetArrayF90','PUMSeqAIJGetArray'])
routine_list.append(['MatSeqAIJRestoreArrayF90','PUMSeqAIJRestoreArray'])
routine_list.append(['MatGetRowIJF90','PUMGetRowIJ'])
routine_list.append(['MatRestoreRowIJF90','PUMRestoreRowIJ'])
routine_list.append(['VecGetArrayF90','PUVGetArray'])
routine_list.append(['VecRestoreArrayF90','PUVRestoreArray'])
routine_list.append(['VecGetArrayReadF90','PUVGetArrayRead'])
routine_list.append(['VecRestoreArrayReadF90','PUVRestoreArrayRead'])
routine_list.append(['VecDuplicateVecsF90','PUVDuplicateVecs'])
routine_list.append(['VecDestroyVecsF90','PUVDestroyVecs'])
routine_list.append(['ISGetIndicesF90','PUIGetIndices'])
routine_list.append(['ISRestoreIndicesF90','PUIRestoreIndices'])

def compile_regex(list_):
    new_list = []
    for entry in list_:
        new_entry = []
        new_entry.append(re.compile(r'\s*call\s*'+str(entry[0])+r'\s*\('))
        for i in range(len(entry)):
            new_entry.append(entry[i])
        new_list.append(new_entry)
    return new_list

def refactor_routine(routine_info,line,fin,fout,debug_info):
    old_name = routine_info[1]
    new_name = routine_info[2]
    regex_string = r'(?P<pre>\s*call\s*)'+str(old_name)+r'\s*\('+\
                   r'(?P<post>.*)'
    regex_ = re.compile(regex_string,flags=myflags)
    m = re.search(regex_,line)
    if not m:
        print(f'error in {debug_info[0]} at line {debug_info[1]} for routine {debug_info[2]}')
        sys.exit(0)
    new_line = m.group('pre').rstrip()+' '+new_name+'('+m.group('post').rstrip()
    fout.write(new_line+'\n')

def refactor(filename,regex_list,overwrite):
    fin = open(filename,'r')
    fout = open(filename+'.tmp','w')
    prev_line = ''
    iline = 0
    for line in fin:
        iline += 1
        found = False
        if line.lstrip().startswith('!'):
            fout.write(line)
            continue
        for routine in regex_list:
            m = re.search(routine[0],line)
            if m:
                refactor_routine(routine,line,fin,fout,
                                 [filename,iline,routine[0]])
                found = True
                break
        if not found:
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
    regex_list = compile_regex(routine_list)
    if test:
        filename_list = ['test']
    else:
        filename_list = get_source_file_roots()
        filename_list.remove('petsc_utility')
    print(filename_list)
    for filename in filename_list:
        filename += '.F90'
        print(f'refactoring {filename}')
        refactor(filename,regex_list,overwrite)


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
