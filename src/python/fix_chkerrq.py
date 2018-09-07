# swap_comment_character.py
import sys
import shutil
import os
import fnmatch
import re

debug_ = False
  
def strip_trailing_comment(line,mo):
  newline = line
  i = newline.find('!')
  if i > -1:
    newline = newline[0:i-1]
  return newline

def parse_file(filename):
  f = open(filename,'r')
  lines = f.readlines()
  f.close

  line_index=0
  for line in lines:
    if line_index>0:
        prev_index=line_index-1
        line_prev=lines[prev_index].lower().strip()
        linelc = line.lower().strip()
        if linelc.find("chkerrq")>0 and line_prev.find("&")>0:
            while line_prev.find("&")>0:
                prev_index-=1
                line_prev=lines[prev_index].lower().strip()
            prev_index+=1
            line_prev=lines[prev_index].lower().strip()
            if line_prev.find("if")==0:
                print("Error in file: "+filename)
                print("Conditional Line "+str(prev_index)+": "+line_prev)
                print("CHKERRQ Line "+str(line_index)+": "+linelc)
    line_index+=1

def get_filenames():
  # Obtain list of source files
  source_file_list = []
  for line in open('pflotran_object_files.txt','r'):
    # find .o file
    # could use re.split() here, but too complicated.
    w = line.split('}')
    if len(w) == 2:
      w2 = w[1].split('.o')
      source_file_list.append(w2[0]+'.F90')
  source_file_list.append('pflotran.F90')
  source_file_list.sort()
  return source_file_list


filenames = get_filenames()
for filename in filenames:
  parse_file(filename)
  

print('done')
