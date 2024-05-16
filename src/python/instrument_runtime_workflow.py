import sys
import shutil
import os
import re

routine_count = 0
routine_tuple = tuple(['abc','def'])
header_tuple = tuple(['implicit none','#','save','end ',
                      'parameter ','parameter(','import ','data '])

end_routine = re.compile(r'end\s*(subroutine|function)')
routine_name = re.compile(r'\s*(subroutine|function)\s*(?P<routine_name>\w+)',re.IGNORECASE)
petsc_statement = re.compile(r'\s*Petsc[a-zA-Z]')
use_statement = re.compile(r'\s*use\s+[a-zA-Z]')
include_statement = re.compile(r'\s*include\s+[a-zA-Z]')
external_statement = re.compile(r'\s*external\s+[a-zA-Z]')
type_statement = re.compile(r'\s*(type|class)\s*\(')
double_colon = re.compile(r'\s*::\s*([a-zA-Z]|&)')
interface_statement = re.compile(r'\s*(interface|abstract\s*interface)')
end_interface_statement = re.compile(r'\s*end\s*interface')
begin_function = re.compile(r'\s*(recursive|pure|petscint|petscreal|petscbool|integer|real\*8|logical)\s*function\s+')
begin_subroutine = re.compile(r'(\s*subroutine|\s*recursive\s+subroutine)\s+')
separate_function = re.compile(r'\s+function\s+')
begin_module = re.compile(r'\s*module\s+')
use_pflotran_constants_statement = re.compile(r'\s*use \s*pflotran_constants_module',re.IGNORECASE)
return_statement = re.compile(r'\.*return\s*$',re.IGNORECASE)

skip_routines_tuple = tuple(['PushCallStack','PopCallStack','PrintCallStack'])
skip_file_roots_tuple = tuple(['pflotran_provenance'])

verbosity_level = 0

def is_module_start(line):
    line2 = line.strip().lower()
    found = False
    if begin_module.match(line2):
        found = True
    return found

def add_missing_use_statement(fin,fout,line,filename):
    if filename.startswith('pflotran_constants'):
        return
    while True:
        line = fin.readline()
        if verbosity_level > 9:
            print(line)
        if use_pflotran_constants_statement.match(line):
            fout.write(line)
            return
        if line.strip().lower().startswith('implicit none'):
            fout.write('use PFLOTRAN_Constants_module\n')
            fout.write(line)
            return
        if line.strip().lower().startswith('contains'):
            print('"contains" found without implicit none in {}'.format(filename))
            sys.exit(1)
        fout.write(line)

def line_ends_with_return(line):
    w = line.split('!') # remove any training comments
    line = w[0]
    return return_statement.search(line)

def insert_pop(line):
    w = line.split('!') # remove any training comments
    line = w[0]
    line2 = line.strip().lower()
    new_line = 'abc'
    if line2.startswith('return') and line2.endswith('return'):
        new_line = 'call PopCallStack();return\n'
    else:
        w = line2.strip().rsplit('return')
        if w[0].strip().endswith(')'):
            new_line = w[0]+' then;call PopCallStack();return;endif;\n'
        else:
            print(line)
            sys.exit(1)
    return new_line

def is_routine_start(line):
    line2 = line.strip().lower()
    found = False
    if begin_subroutine.match(line2):
        found = True
    elif line2.startswith('function'):
        found = True
    elif 'function' in line2:
        if begin_function.match(line2):
            found = True
            if line2.startswith('pure'):
                found = False
        elif not separate_function.search(line2) == None:
            if verbosity_level > 0:
                print('Possible function miss: {}'.format(line))
    if found:
        if get_routine_name(line) in skip_routines_tuple:
            found = False
    return found

def first_executable_line(line):
    line_lower = line.strip().lower()
    first_exec_line = True
    if len(line_lower) == 0:
        first_exec_line = False
    elif line_lower.startswith(header_tuple):
        first_exec_line = False
    elif double_colon.search(line):
        first_exec_line = False
    elif petsc_statement.match(line):
        first_exec_line = False
    elif include_statement.match(line_lower):
        first_exec_line = False
    elif use_statement.match(line_lower):
        first_exec_line = False
    elif type_statement.match(line_lower):
        first_exec_line = False
    elif external_statement.match(line_lower):
        first_exec_line = False
    if first_exec_line:
        if verbosity_level > 0:
            print(line)
    return first_exec_line

def line_ends_with_ampersand(line):
    return line.find('&') > -1 and \
           (line.find('!') == -1 or
            line.find('!') > line.find('&'))

def is_routine_end(line):
    line = line.lower()
    return end_routine.match(line.strip()) != None

def get_routine_name(line):
    m = routine_name.search(line)
    if not m:
        print('get_routine_name issue: {}'.format(line))
        sys.exit(1)
    name = m.group('routine_name')
    return name

def check_preprocessing_conditional(fin,fout,line):
    if line.strip().lower().startswith('#if'):
        fout.write(line)
        while True:
            line = fin.readline()
            if check_preprocessing_conditional(fin,fout,line):
                continue
            fout.write(line)
            if line.strip().lower().startswith('#endif'):
                break
        return True
    else:
        return False

def refactor_routine(fin,fout,line):
    global routine_count
    routine_count += 1
    routine_name = get_routine_name(line)
    prev_line_end_with_ampersand = line_ends_with_ampersand(line)
    if verbosity_level > 0:
        print(routine_name)
    flag = False
    while True:
        line = fin.readline()
        if verbosity_level > 9:
            print(line)
        if check_preprocessing_conditional(fin,fout,line):
            continue
        if interface_statement.match(line):
            fout.write(line)
            skip_interface(fin,fout)
            continue
        if len(line) == 0:
            print('Premature end of routine: {}'.format(line))
            sys.exit(1)
        if line.strip().startswith('!'):
            fout.write(line)
            continue
        if not flag and first_executable_line(line) and \
           not prev_line_end_with_ampersand:
            flag = True
            fout.write('call PushCallStack("{}")\n'.format(routine_name))
        routine_end = is_routine_end(line)
        if routine_end:
            fout.write('call PopCallStack()\n')
        if line_ends_with_return(line):
            fout.write(insert_pop(line))
        else:
            fout.write(line)
        if not flag:
            prev_line_end_with_ampersand = line_ends_with_ampersand(line)
        if routine_end:
            return

def skip_interface(fin,fout):
    while True:
        line = fin.readline()
        if verbosity_level > 9:
            print(line)
        fout.write(line)
        if end_interface_statement.match(line):
            break

def instrument_file(filename):
    fin = open(filename,'r')
    fout = open(filename+'.tmp','w')
    while True:
        line = fin.readline()
        if verbosity_level > 9:
            print(line)
        if len(line) == 0:
            break
        fout.write(line)
        if interface_statement.match(line):
            skip_interface(fin,fout)
            continue
        if is_module_start(line):
            add_missing_use_statement(fin,fout,line,filename)
            continue
        if is_routine_start(line):
            refactor_routine(fin,fout,line)
    fin.close()
    fout.close()

instrument_headers = '''

  character(len=MAXWORDLENGTH) :: call_stack(1000)
  PetscInt :: call_stack_count = 0
  public :: PushCallStack, &
            PopCallStack

'''

instrument_routines = '''
! ************************************************************************** !

subroutine PushCallStack(routine_name)
  !
  ! Appends routine name to a long string
  !
  ! Author: Glenn Hammond
  ! Date: 04/15/24
  !
  implicit none

  character(len=*) :: routine_name

  character(len=:), allocatable :: string
  PetscInt :: i

  call_stack_count = call_stack_count + 1
  call_stack(call_stack_count) = routine_name

  string = ''
  do i = 1, call_stack_count
    if (i > 1) then
      string = string // ','
    endif
    string = string // trim(call_stack(i))
  enddo
  print *, 'Call Stack: ' // string

end subroutine PushCallStack

! ************************************************************************** !

subroutine PopCallStack()
  !
  ! Appends routine name to a long string
  !
  ! Author: Glenn Hammond
  ! Date: 04/15/24
  !
  implicit none

  call_stack_count = call_stack_count - 1
  if (call_stack_count < 0) then
    print *, 'call stack count < 0 ', call_stack_count
  endif

end subroutine PopCallStack

! ************************************************************************** !

subroutine PrintCallStack()
  !
  ! Appends routine name to a long string
  !
  ! Author: Glenn Hammond
  ! Date: 04/15/24
  !
  implicit none

  character(len=:), allocatable :: string
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: i


  string = ''
  do i = 1, call_stack_count
    if (i > 1) then
      string = string // ','
    endif
    string = string // trim(call_stack(i))
  enddo
  write(word,*) call_stack_count
  print *, 'Call Stack [' // trim(adjustl(word)) // ']: ' // string

end subroutine PrintCallStack

'''

def instrument_pflotran_constants():
    filename = 'pflotran_constants.F90'
    fin = open(filename,'r')
    fout = open(filename+'.tmp','w')
    while True:
        line = fin.readline()
        if len(line) == 0:
            break
        if line.strip().startswith('contains'):
            fout.write(instrument_headers)
        if line.strip().startswith('end module PFLOTRAN_Constants_module'):
            fout.write(instrument_routines)
        fout.write(line)
    fin.close()
    fout.close()
    shutil.copy(filename+'.tmp',filename)

def get_filename(root,suffix):
    filename = []
    filename.append(root)
    filename.append(suffix)
    filename = '.'.join(filename)
    return filename

#def instrument_pflotran_constants():
#    filename = 'pflotran_constants.F90'
#    fin = open(filename,'r')
#    fout = open(filename+'.tmp','w')
#    while True:
#        line = fin.readline()
#        if len(line) == 0:
#            break
#        if line.strip().startswith('end module PFLOTRAN_Constants_module'):
#
#        fout.write(line)

#filename = 'temp.F90'
#instrument_file(filename)

# Obtain list of source files
source_file_roots = []
for line in open('pflotran_object_files.txt','r'):
    # find .o file
    # could use re.split() here, but too complicated.
    w = line.split('}')
    if len(w) == 2:
        w2 = w[1].split('.o')
        source_file_roots.append(w2[0])
#source_file_roots = ['pflotran_constants']
source_file_roots.append('pflotran')

# Alphabetize
source_file_roots.sort()
if verbosity_level > 0:
    print(source_file_roots)
#os.remove('pflotran_provenance.F90')
file_count = 0
for root in source_file_roots:
    if root in skip_file_roots_tuple:
        continue
    file_count += 1
    if not root.endswith('pf'):
      filename = get_filename(root,'F90')
    else:
      filename = root

 #   print(filename)
    instrument_file(filename)
    shutil.copy(filename+'.tmp',filename)

instrument_pflotran_constants()

os.system('python3 ../python/pflotran_dependencies.py')

print(file_count)
print('done!')


