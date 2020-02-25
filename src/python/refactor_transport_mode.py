# swap_comment_character.py
import sys
import shutil
import os
import fnmatch

def refactor_options(f,f2):
    options_block_started = False
    skip_file_flag = False
    while True:
        line = f.readline()
        string = line.strip().upper()
        if string.startswith('MODE'):
            # return as file has already been fixed.
            # flag causes the script to skip file
            skip_file_flag = True
            break
        elif string.startswith('/') or string.startswith('END'):
            if options_block_started:
                new_line = '      /\n'
                f2.write(new_line)
            f2.write(line)
            break
        elif string.startswith('GLOBAL_IMPLICIT'):
            continue
        elif string.startswith('NUMERICAL_JACOBIAN'):
            new_line = '      OPTIONS\n'
            f2.write(new_line)
            new_line = '        NUMERICAL_JACOBIAN\n'
            f2.write(new_line)
            new_line = '      /\n'
            f2.write(new_line)
        else:
            if not options_block_started:
                options_block_started = True
                new_line = '      OPTIONS\n'
                f2.write(new_line)
            f2.write('  '+line.rstrip()+'\n')
    return skip_file_flag


def refactor_file(filename):
    f = open(filename,'r')
    f2 = open(filename+'.tmp','w')
    flag = False
    skip_file_flag = False
    while True:
        if skip_file_flag:
            break
        line = f.readline()
        if not line:
            break
        string = line.strip().upper()
        if string.startswith('SUBSURFACE_TRANSPORT'):
            flag = True
            f2.write(line)
            new_line = '      MODE GIRT\n'
            f2.write(new_line)
            skip_file_flag = refactor_options(f,f2)
            continue
        if string.startswith('NUCLEAR_WASTE_TRANSPORT'):
            flag = True
            f2.write(line.replace('NUCLEAR_WASTE_TRANSPORT','SUBSURFACE_TRANSPORT'))
            new_line = '      MODE NWT\n'
            f2.write(new_line)
            skip_file_flag = refactor_options(f,f2)
            continue
        f2.write(line)
    f.close()
    f2.close()
    if flag and not skip_file_flag:
        print('File {} has been updated.'.format(filename))
        # using shutil.move adds ^M to end of lines.
        os.remove(filename)
        shutil.copy(filename+'.tmp',filename)
    os.remove(filename+'.tmp')

suffix = '*.in'
for root, dirnames, filenames in os.walk('.'):
    for filename in fnmatch.filter(filenames,suffix):
        filename = os.path.join(root,filename)
#        print(filename)
        refactor_file(filename)

print('done')
