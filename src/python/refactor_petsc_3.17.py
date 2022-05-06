import sys
import shutil
import os
import re

petsc_prefix_list = ['Petsc', 
                     'Vec',
                     'Mat',
                     'DM',
                     'AO',
                     'IS',
                     'TS',
                     'SNES',
                     'KSP',
                     'PC',
                     'MPI']

not_petsc_prefix_list = ['Material',
                         'Matrix',
                         'TSPre',
                         'TSrc',
                         'TS_',
                         'MatGetSub',
                         'Isoth',
                         'pck']

path = ''

test = True
test = False
overwrite_files = True

petsc_call = False
#petsc_call = True

error_check_string = 'PetscCall'
error_check_string = 'CHKERRQ'

enforce_ierr = False

#overwrite_files = False

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

myflags = re.M|re.I

string = '\s*call\s*(?:'+'|'.join(petsc_prefix_list)+')'
not_string = '\s*call\s*(?:'+'|'.join(not_petsc_prefix_list)+')'
#petsc_function = re.compile(r'call\s*(?:Mat|Vec)',flags=myflags)
petsc_function = re.compile(string,flags=myflags)
non_petsc_function = re.compile(not_string,flags=myflags)

#petsc_swap = re.compile(r'(?P<indentation>\s*)'
#                        r'call\s*(?P<function_name>\w+)?\s*'
#                        r'[(]\.*[)]',
#                        flags=myflags)

petsc_swap = re.compile(r'(?P<indentation>\s*)'
                        r'call\s*(?P<function_name>\w+)?\s*'
                        r'[(]\s*(?P<arguments>.*)[)]\s*',
                        flags=re.I)

chkerrq_swap = re.compile(r'\s*CHKERRQ\s*[(]\s*(?P<arguments>.*)\s*[)]\s*',
                          flags=re.I)



#petsc_swap = re.compile(r'Vec',
#                        flags=myflags)


def is_petsc_function(line):
    if petsc_function.match(line):
        if not non_petsc_function.match(line):
            return True
    return False

def get_petsc_error_check(arg):
    line = ';'
    line += error_check_string
    line += '('+arg.strip()+')'
    return line

def chkerrq_sub(match):
    arg = match.group('arguments')
    if not arg.startswith('ierr'):
        print('CHKERRQ arg not ierr')
    if enforce_ierr:
        arg = 'ierr'
    return get_petsc_error_check(arg)

def function_sub(match):
    indentation = len(match.group('indentation'))
    name = match.group('function_name').strip()
    arguments = match.group('arguments').split(',')
    list = []
    list.append(' '*indentation)
    if petsc_call:
        list.append('PetscCall(')
    else:
        list.append('call ')
    list.append(name)
    list.append('(')
    args = []
    for arg in arguments:
        args.append(arg.replace(' ',''))
    arguments = ','.join(args)
    list.append(arguments)
    list.append(')')
    if petsc_call:
        list.append(')')
    return ''.join(list)

def count_parentheses(string):
    icount = 0
    s = string.split('!')
    if len(s) > 1:
        print('Comments in subroutine argument list')
    string = s[0]
    string = string.rstrip(' &')
    for char in string:
        if char == '(':
            icount += 1
        elif char == ')':
            icount -= 1
    return icount, string
    
def indentation_parentheses(string):
    string = string.split('!')[0]
    icount = 0
    parentheses_count = 0
    for char in string:
        if char == '(':
            parentheses_count += 1
            if petsc_call:
                if parentheses_count == 2:
                    break
            else:
                break
        icount += 1
    icount += 1
    return icount

def refactor_petsc_function(f,line):
#    print(line)
#    line = petsc_swap.sub(function_sub,line)
#    line = petsc_function.sub(function_sub,line)
    missing_chkerr = False
    line = line.rstrip()
    print(line)
    i,s = count_parentheses(line)
    while i > 0:
        line = f.readline().rstrip()
        print(line)
        ii,ss = count_parentheses(line)
        i += ii
        s += ss.strip()
        if i < 0:
            sys.exit('Error counting parentheses')
    s = s.split(';')
    line = petsc_swap.sub(function_sub,s[0])
    print(line)
    if not petsc_call:
        if len(s) > 1:
            line += chkerrq_swap.sub(chkerrq_sub,s[1])
        else:
            missing_chkerr = True
            args = line.split(',')
            if len(args) > 1:
                arg = args[-1].split(')')[0]
            else:
                arg = line.split('(')[-1].split(')')[0]
            line += get_petsc_error_check(arg)
    indentation = ' '*indentation_parentheses(line)
    indent_count = len(indentation)
    print(line)
    # break line into call statement and arguments
    w0 = line.split(',')
    w1 = w0[0].split('(',1)
    w1[0] = w1[0].rstrip()+'('
    w = w1+w0[1:]
    # if any of the words have an open or close parenthese without a
    # match, we have to merge
    for i in range(1,len(w)):
        icount, s = count_parentheses(w[i])
        j = i+1
        if icount > 0 and j < len(w):
            jcount = 0
            new_word = w[i]
            while not icount + jcount == 0:
                jcount, s = count_parentheses(w[j])
                new_word += ',' + w[j]
                w[j] = ''
            w[i] = new_word
    # remove empty words
    ww = []
    for word in w:
        if len(word) > 0:
            ww.append(word)
    w = ww
    line = w[0]
    icount = len(w[0])
    max_char = 76
    wrap_flag = True
    for word in w[1:]:
        lenw = len(word)
        if icount + lenw > max_char:
            max_char = 76
            if not wrap_flag:
                line += ', &\n'+ indentation
            icount = indent_count
            # still beyond
            if icount + lenw > max_char:
                w = word.split('%')
                if len(w) > 1:
                    overage = icount + lenw - max_char
                    i = len(w)-1
                    while(overage > 0 and i > -1):
                        overage -= len(w[i])
                        i -= 1
                    word = '%'.join(w[0:i+1])
                    line += word + '% &\n'+ indentation
                    line += '  ' + '%'.join(w[i+1:]) + ', &\n'+ indentation
                    icount = indent_count
                    wrap_flag = True
                    continue
                else:
                    deindent = 0
                    for i in range(2,8,2):
                        if icount + lenw <= max_char + i:
                            deindent = i
                    if icount + lenw > max_char + deindent:
                        ichar = icount + lenw - max_char
                        ichar -= indent_count
                        word1 = word[:ichar]
                        word2 = word[ichar:]
                        line += word1 + '&\n' + indentation
                        line += '  &' + word2
                        icount = indent_count + len(word2)
                    else:
                        line = line[:-deindent]
                        line += word + ', &\n'+ indentation
                        wrap_flag = True
                    continue
        else:
            if not wrap_flag:
                line += ','
                icount += 1
            wrap_flag = False
        line += word
        icount += lenw
#        print(line)
    return line,missing_chkerr


    
#line = '  call   VecZeroEntries ( res_pert, ierr) ; CHKERRQ(ierr) ' 
#line = '  call   VecZeroEntries ( res_pert, ierr)' 
#print(is_petsc_function(line)) 
'''
filename = 'test.txt'
f = open_f(filename)
line = f.readline()
new_line = refactor_petsc_function(f,line)
f.close()
print(new_line)

sys.exit(0)
'''

def refactor_file(filename):
    num_petsc_functions = 0
    num_unchanged = 0
    num_updates = 0
    num_missing_chkerrq = 0
    f = open_f(filename)
    if overwrite_files:
        f2 = open_f2(filename)
#    content = f.read()
#    revised_content = input_type.sub('class(input_type)',content)
#    f2.write(revised_content)  
    while(True):
        line = f.readline()
        if len(line) < 1:
            break
        if is_petsc_function(line):
            num_petsc_functions += 1
            orig_line = line.rstrip('\n')
            line,missing_chkerr = refactor_petsc_function(f,line)
            if line == orig_line:
                num_unchanged += 1
            else:
                num_updates += 1
            if missing_chkerr:
                num_missing_chkerrq += 1
            line = line.rstrip()+'\n'
            print(line)
        if overwrite_files:
            f2.write(line)
    f.close()
    if overwrite_files:
        f2.close()
        if not test:
            mv_file(filename)
    return num_petsc_functions, num_unchanged, num_updates, num_missing_chkerrq

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

# Alphabetize
source_file_roots.sort()
print(source_file_roots)
#os.remove('pflotran_provenance.F90')
if test:
    source_file_roots=['test']
file_count = 0
num_petsc_functions = 0
num_unchanged = 0
num_updates = 0
num_missing_chkerr = 0
for root in source_file_roots:
    file_count += 1
    if not root.endswith('pf'):
        filename = get_filename(root,'F90')
    else:
        filename = root
#    print(filename)
    n0, n1, n2, n3 = refactor_file(filename)
    num_petsc_functions += n0
    num_unchanged += n1
    num_updates += n2
    num_missing_chkerr += n3

print('Number of files: {}'.format(file_count))
print('Number of PETSc functions: {}'.format(num_petsc_functions))
print('Number unchanged {}'.format(num_unchanged))
print('Number of updates {}'.format(num_updates))
print('Number of missing CHKERRQ: {}'.format(num_missing_chkerr))
print('done!')
