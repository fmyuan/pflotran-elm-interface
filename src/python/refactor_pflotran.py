import sys
import shutil
import os
import re

test = True
test = False
overwrite_files = True
if test:
  overwrite_files = False
move_to_build = True
move_to_build = False

use_option_class_procedures = False
use_input_class_procedures = True

add_class_procedures = True

use_option_class = False
use_input_class = True

debug_printmsgsub = True
if not test:
  debug_printmsgsub = False
debug_inputsub = True
if not test:
  debug_inputsub = False
if move_to_build:
  path = 'build/'
else:
  path = './'
  
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
  if overwrite_files:
    # using shutil.move adds ^M to end of lines.
    shutil.copy(path+filename+'.tmp',path+filename)
    os.remove(path+filename+'.tmp')

def cleanup(line):
  return line.rstrip()
  # the below doesnt work 
  w = line.split('!')
  comment = ''
  if len(w) == 2:
    comment = ' !' + w[1]
  base = w[0]
  # remove spaces beween commas
  base = re.sub(r'\s*,\s',',',base)
  if base.endswith('&'):
    # insert space before &
    base = re.sub(r'(.*)\&(.*)',r'\1 &\2',base)
  return base + comment 
  
myflags = re.M|re.I
# input
input_search = re.compile(r'input',flags=myflags)
input_type = re.compile(r'type\s*[(]\s*input_type\s*[)]',flags=myflags)

error = re.compile(r'\s*inputerror\s*[(]\s*'
                   r'(?P<input_object>input\w*)\s*',
                   flags=myflags)
                   
input_ = re.compile(r'call\s*(?P<routine>input\w*)\s*[(]\s*'
                      r'(?P<input_object>.*?input\w*)' # have to have input
                      r'(?P<input_member>[%]\w+)?'
                      r'(?P<comma>\s*,\s*)?'
                      r'(?P<arguments>.*)?',
                      flags=myflags)
                
# option            
option_search = re.compile(r'option',flags=myflags)
#insert_option = re.compile(r'subroutine\s*print',flags=myflags)
insert_option = re.compile(r'\s{1}p(?P<routine>rint\w*msg\w*)',flags=myflags)

option_type = re.compile(r'type\s*[(]\s*option_type\s*[)]',flags=myflags)

#                      r'(?P<option_object>(\w+[%]?)+)?'
printmsg = re.compile(r'call\s*p(?P<routine>rint\w*msg\w*)\s*[(]\s*'
                      r'(?P<option_object>.*?option)' # have to have option
                      r'(?P<option_member>[%]\w+)?'
                      r'(?P<comma>\s*,\s*)?'
                      r'(?P<arguments>.*)?',
                      flags=myflags)
printmsgf = re.compile(r'call\s*print\w*msg\w*\s*[(]\s*'
                       r'(?P<option_object>\w+)?\s*\W\s*',
                       flags=myflags)

def printmsgsub(match):
  option_object = match.group('option_object')
  routine = match.group('routine')
  option_member = match.group('option_member')
  comma = match.group('comma')
  arguments = match.group('arguments')
  if debug_printmsgsub:
    print('printmsgsub')
    print('routine:   '+routine+'EOL')
    if option_object:
      print('object:    '+option_object+'EOL')
    else:
      print('object:    nada'+'EOS')
    if option_member:
      print('member:    '+option_member+'EOL')
    else:
      print('member:    nada'+'EOS')
    if comma:
      print('comma:     '+comma+'EOL')
    else:
      print('comma:     nada'+'EOL')
    if arguments:
      print('arguments: '+arguments+'EOL')
    else:
      print('arguments: nada'+'EOL')
  # remove Input from the routine name
  routine = re.sub(r'^option','',routine,flags=myflags)
  if routine and option_object and arguments:
    arguments = cleanup(arguments)
    if option_member:
      if debug_printmsgsub:
        print('case1')
      return 'call {}%P{}(option{}{}'.format(option_object, \
                                             routine, \
                                             option_member, \
                                             arguments)
    else:
      if debug_printmsgsub:
        print('case2')
      return 'call {}%P{}({}'.format(option_object, \
                                     routine, \
                                     arguments)
  elif routine and option_object:
    if option_member:
      if debug_printmsgsub:
        print('case3')
      return 'call {}%P{}(option{}'.format(option_object, \
                                           routine, \
                                           option_member)
    else:
      if debug_printmsgsub:
        print('case4')
      return 'call {}%P{}('.format(option_object, \
                                   routine)
  elif routine:
    if debug_printmsgsub:
      print('case5')
    return 'call option%P{}('.format(routine)
  return 'bogus'
  
def inputsub(match):
  input_object = match.group('input_object')
  routine = match.group('routine')
  input_member = match.group('input_member')
  comma = match.group('comma')
  arguments = match.group('arguments')
  if debug_inputsub:
    print('inputsub')
    print('routine:   '+routine+'EOL')
    if input_object:
      print('object:    '+input_object+'EOL')
    else:
      print('object:    nada'+'EOS')
    if input_member:
      print('member:    '+input_member+'EOL')
    else:
      print('member:    nada'+'EOS')
    if comma:
      print('comma:     '+comma+'EOL')
    else:
      print('comma:     nada'+'EOL')
    if arguments:
      print('arguments: '+arguments+'EOL')
    else:
      print('arguments: nada'+'EOL')
  reconstruct = False
  if input_object:
    if not input_object.startswith('input'):
      reconstruct = True
  # remove Input from the routine name
  routine = re.sub(r'^input','',routine,flags=myflags)
  if not routine.startswith(('Read',
                             'Error',
                             'Default',
                             'CheckExit',
                             'CheckMandatoryUnits',
                             'FindStringErrorMsg')):
    reconstruct = True
  elif routine.startswith(('ReadPflotran',
                           'ReadToBuffer')):
    reconstruct = True
  if reconstruct:
    line = 'call Input{}('.format(routine)
    if input_object:
      line += '{}'.format(input_object)
    if input_member:
      line += '{}'.format(input_member)
    if comma:
      line += ','
    if arguments:
      line += '{}'.format(arguments)
    return line
  if routine and input_object and arguments:
    arguments = cleanup(arguments)
    if input_member:
      if debug_inputsub:
        print('case1')
      return 'call {}%{}({}{},{}'.format(input_object, \
                                        routine, \
                                        input_object, \
                                        input_member, \
                                        arguments)
    else:
      if debug_inputsub:
        print('case2')
      return 'call {}%{}({}'.format(input_object, \
                                    routine, \
                                    arguments)
  elif routine and input_object:
    if input_member:
      if debug_inputsub:
        print('case3')
      return 'call {}%{}{}({}'.format(input_object, \
                                      routine, \
                                      input_object, \
                                      input_member)
    else:
      if debug_inputsub:
        print('case4')
      return 'call {}%{}('.format(input_object, \
                                  routine)
  elif routine:
    if debug_inputsub:
      print('case5')
    return 'call input%{}('.format(routine)
  return 'bogus'

#print(re.sub('.*option','replaced','ctx%option'))
#sys.exit(1)
#m = re.match('\.*option','ctx%option')
#print(m.groups)

#line3 = 'call printErrMsg(ctx%option,ctx%option%io_buffer)'
#print(line3)
#print(printmsg.sub(printmsgsub,line3))
#sys.exit(1)          

#input_type = re.compile(r'type\s*[(]\s*input_type',re.IGNORECASE)

def revise_option():
  print('refactoring option')
  filename = 'option.F90'
  f = open_f(filename)
  f2 = open_f2(filename)
  for line in f:
    if line.strip().startswith('end type option_type') and \
       add_class_procedures:
#    procedure, nopass :: PrintMsg2 => OptionPrintMsg2
#    procedure, nopass :: PrintMsgAnyRank2 => OptionPrintMsgAnyRank2
#    procedure, nopass :: PrintMsgByRank2 => OptionPrintMsgByRank2
#    procedure, nopass :: PrintErrMsgByRank2 => OptionPrintErrMsgByRank2
#    procedure, nopass :: PrintErrMsgNoStopByRank2 => &
#                           OptionPrintErrMsgNoStopByRank2
#    procedure, nopass :: PrintErrMsg2 => OptionPrintErrMsg2
#    procedure, nopass :: PrintWrnMsg2 => OptionPrintWrnMsg2
      f2.write('''  contains
    procedure :: PrintMsg1 => OptionPrintMsg1
    procedure :: PrintMsg2 => OptionPrintMsg2
    generic, public :: PrintMsg => PrintMsg1,PrintMsg2
    procedure :: PrintMsgAnyRank1 => OptionPrintMsgAnyRank1
    procedure, nopass :: PrintMsgAnyRank2 => OptionPrintMsgAnyRank2
    generic, public :: PrintMsgAnyRank => PrintMsgAnyRank1, &
                                          PrintMsgAnyRank2
    procedure :: PrintMsgByRank1 => OptionPrintMsgByRank1
    procedure :: PrintMsgByRank2 => OptionPrintMsgByRank2
    generic, public :: PrintMsgByRank => PrintMsgByRank1, &
                                         PrintMsgByRank2
    procedure :: PrintErrMsgByRank1 => OptionPrintErrMsgByRank1
    procedure :: PrintErrMsgByRank2 => OptionPrintErrMsgByRank2
    generic, public :: PrintErrMsgByRank => PrintErrMsgByRank1, &
                                            PrintErrMsgByRank2
    procedure :: PrintErrMsgNoStopByRank1 => OptionPrintErrMsgNoStopByRank1
    procedure :: PrintErrMsgNoStopByRank2 => &
                           OptionPrintErrMsgNoStopByRank2
    generic, public :: PrintErrMsgNoStopByRank => &
                           PrintErrMsgNoStopByRank1, &
                           PrintErrMsgNoStopByRank2
    procedure :: PrintErrMsg1 => OptionPrintErrMsg1
    procedure :: PrintErrMsg2 => OptionPrintErrMsg2
    generic, public :: PrintErrMsg => PrintErrMsg1,PrintErrMsg2
    procedure :: PrintWrnMsg1 => OptionPrintWrnMsg1
    procedure :: PrintWrnMsg2 => OptionPrintWrnMsg2
    generic, public :: PrintWrnMsg => PrintWrnMsg1,PrintWrnMsg2
    procedure :: PrintErrMsgToDev => OptionPrintErrMsgToDev
    procedure :: PrintErrMsgByRankToDev => OptionPrintErrMsgByRankToDev
    procedure :: PrintMsgByCell => OptionPrintMsgByCell
    procedure :: PrintVerboseMsg => OptionPrintVerboseMsg
    procedure :: InitMPI1 => OptionInitMPI1
    procedure :: InitMPI2 => OptionInitMPI2
    generic, public :: InitMPI => InitMPI1,InitMPI2\n''')
#    line = insert_option.sub('subroutine OptionPrint',line)
    if use_option_class_procedures or add_class_procedures:
      line = insert_option.sub(' OptionP\g<routine>',line)
    f2.write(line)
  f.close()
  f2.close()
  mv_file(filename)
    

def revise_input_aux():
  print('refactoring input_aux')
  filename = 'input_aux.F90'
  f = open_f(filename)
  f2 = open_f2(filename)
  for line in f:
    if line.strip().startswith('end type input_type') and \
       add_class_procedures:
      f2.write('''  contains
    procedure :: ReadWord1 => InputReadWord1
    procedure, nopass :: ReadWord2 => InputReadWord2
    generic, public :: ReadWord => ReadWord1,ReadWord2
    
    
    
    procedure :: ReadNChars1 => InputReadNChars1
    procedure, nopass :: ReadNChars2 => InputReadNChars2
    generic, public :: ReadNChars => ReadNChars1,ReadNChars2

    procedure, public :: ReadWordDbaseCompatible => &
                           InputReadWordDbaseCompatible
    procedure, public :: ReadQuotedWord => InputReadQuotedWord
    
    procedure :: ReadInt1 => InputReadInt1
    procedure, nopass :: ReadInt2 => InputReadInt2
#if defined(PETSC_USE_64BIT_INDICES) && (PETSC_SIZEOF_MPI_FINT * PETSC_BITS_PER_BYTE != 64)
    procedure :: ReadInt3 => InputReadInt3
    procedure, nopass :: ReadInt4 => InputReadInt4
    generic, public :: ReadInt => ReadInt1,ReadInt2,ReadInt3,ReadInt4
#else
    generic, public :: ReadInt => ReadInt1,ReadInt2
#endif
    
    procedure :: ReadDouble1 => InputReadDouble1
    procedure, nopass :: ReadDouble2 => InputReadDouble2
    generic, public :: ReadDouble => ReadDouble1,ReadDouble2
    
    procedure :: ReadNDoubles1 => InputReadNDoubles1
    procedure, nopass :: ReadNDoubles2 => InputReadNDoubles2
    generic, public :: ReadNDoubles => ReadNDoubles1,ReadNDoubles2
    
    procedure, public :: ErrorMsg1 => InputErrorMsg1
    procedure, public :: ErrorMsg2 => InputErrorMsg2
    generic, public :: ErrorMsg => ErrorMsg1,ErrorMsg2
    
    procedure, public :: ReadStringErrorMsg1 => InputReadStringErrorMsg1
    procedure, public :: ReadStringErrorMsg2 => InputReadStringErrorMsg2
    generic, public :: ReadStringErrorMsg => ReadStringErrorMsg1, &
                                             ReadStringErrorMsg2
    
    procedure, public :: DefaultMsg1 => InputDefaultMsg1
    procedure, public :: DefaultMsg2 => InputDefaultMsg2
    generic, public :: DefaultMsg => DefaultMsg1,DefaultMsg2
    
    procedure :: Error1 => InputError1
    procedure, nopass :: Error2 => InputError2
    generic, public :: Error => Error1,Error2
    
    procedure, public :: FindStringErrorMsg => InputFindStringErrorMsg
    
    procedure, public :: ReadFilename => InputReadFilename   
    
    procedure, public :: CheckExit => InputCheckExit
    
    procedure, public :: ReadAndConvertUnits => InputReadAndConvertUnits
    procedure, public :: CheckMandatoryUnits => InputCheckMandatoryUnits
    
\n''')
    f2.write(line)
  f.close()
  f2.close()
  mv_file(filename)

def refactor_classes(filename):
  f = open_f(filename)
  f2 = open_f2(filename)
#  content = f.read()
#  revised_content = input_type.sub('class(input_type)',content)
#  f2.write(revised_content)  
  for line in f:
    if re.search(option_search,line) and use_option_class:
      line = option_type.sub('class(option_type)',line)
    if re.search(printmsgf,line) and use_option_class:
      if use_option_class_procedures:
        line = printmsg.sub(printmsgsub,line)
      elif add_class_procedures:
        line = insert_option.sub(' OptionP\g<routine>',line)
    if re.search(input_search,line):
      if use_input_class:
        line = input_type.sub('class(input_type)',line)
      if use_input_class_procedures:
        line = error.sub('\g<input_object>%Error(',line)
        line2 = input_.sub(inputsub,line)
        if len(line2) > 1:
          line = line2
    f2.write(line)
  f.close()
  f2.close()
  mv_file(filename)

# Obtain list of source files
source_file_roots = []
for line in open('pflotran_object_files.txt','r'):
  # find .o file
  # could use re.split() here, but too complicated.
  w = line.split('}')
  if len(w) == 2:
    w2 = w[1].split('.o')
#    print(w2[0])
    source_file_roots.append(w2[0])
source_file_roots.append('pflotran')
source_file_roots.append('pflotran_rxn')
source_file_roots.append('pflotran_derivative')
source_file_roots.append('unittests/test_characteristic_curves.pf')
source_file_roots.append('unittests/test_eos_gas.pf')
source_file_roots.append('unittests/test_eos_water.pf')
source_file_roots.append('unittests/test_geometry.pf')
source_file_roots.append('unittests/test_material.pf')
source_file_roots.append('unittests/test_saturation_function.pf')
source_file_roots.append('unittests/test_string.pf')
source_file_roots.append('unittests/test_utility.pf')
source_file_roots.append('../clm-pflotran/clm_pflotran_interface_data')
source_file_roots.append('../clm-pflotran/mapping')
source_file_roots.append('../clm-pflotran/pflotran_interface_main')
source_file_roots.append('../clm-pflotran/pflotran_model')

# Alphabetize
source_file_roots.sort()
print(source_file_roots)
#os.remove('pflotran_provenance.F90')
if test:
  source_file_roots=['test']
if move_to_build:
  copy_files(source_file_roots)
file_count = 0
for root in source_file_roots:
  file_count += 1
  if not root.endswith('pf'):
    filename = get_filename(root,'F90')
  else:
    filename = root
  print(filename)
  refactor_classes(filename)
if use_option_class:
  revise_option()
if use_input_class:
  revise_input_aux()

print(file_count)
print('done!')
