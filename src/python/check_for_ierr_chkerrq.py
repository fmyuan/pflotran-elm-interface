import sys
import shutil
import os
import fnmatch
import re


myflags = re.I

regex_chkerr = re.compile(r'(?P<pre>\w*)\s*\).?CHKERRQ\s*\(\s*(?P<arg>\w*)',flags=myflags)
regex_pre = re.compile(r'(?P<pre>\w*)\s*\).?$')
regex_arg = re.compile(r'^\s*CHKERRQ\s*\(\s*(?P<arg>\w*)')

def check(filename):
    f = open(filename,'r')
    line_count = 0
    prev_line = ''
    for line in f:
        line_count += 1
        m = re.search(regex_chkerr,line)
        if m:
            if not m.group('pre').startswith(m.group('arg')):
                print('Error at line {} of {}:'.format(line_count,filename))
                print('  ({}) : {}'.format(line_count,line.rstrip()))
                print('"{}" versus "{}"\n'.format(m.group('pre'),m.group('arg')))
        else:
            m1 = re.search(regex_pre,prev_line)
            m2 = re.search(regex_arg,line)
            if m1 and m2:
                if not m1.group('pre').startswith(m2.group('arg')):
                    print('Error at line {} of {}:'.format(line_count-1,filename))
                    print('  ({}) : {}'.format(line_count-1,prev_line.rstrip()))
                    print('  ({}) : {}'.format(line_count,line.rstrip()))
                    print('"{}" versus "{}"\n'.format(m1.group('pre'),m2.group('arg')))
        prev_line = line
    f.close()



def main():
    print('\nChecking for missing ierr arguments in PETSc calls.\n')
    suffix = '*.F90'
    for root, dirnames, filenames in os.walk('.'):
        for filename in fnmatch.filter(filenames,suffix):
            filename = os.path.join(root,filename)
#            print(filename)
            check(filename)


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