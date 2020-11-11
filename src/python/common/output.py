#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description: Output class that allows one to print information to screen
             and file simultaneously.
Author: Glenn Hammond
"""
import sys
import os

try:
    pflotran_dir = os.environ['PFLOTRAN_DIR']
except KeyError:
    print('PFLOTRAN_DIR must point to PFLOTRAN installation directory and '+
          'be defined in system environment variables.')
    sys.exit(1)
sys.path.append(pflotran_dir + '/src/python')

from common.error_messaging import *

class Output():
    
    def __init__(self,print_to_screen=True,filename=''):
        self.print_to_screen = print_to_screen
        self.filename = filename
        self.fid = None
        if len(self.filename) > 0:
            self.fid = open(filename,'w')
        self.print_to_screen = print_to_screen
        self.indentation = ''

    def set_indentation(self,num_spaces):
        self.indentation = ' ' * num_spaces
        
    def add_indentation(self,num_spaces):
        self.indentation += ' ' * num_spaces
        
    def subtract_indentation(self,num_spaces):
        l = len(self.indentation)
        self.indentation = ' ' * (l-num_spaces)
        
    def clear_indentation(self):
        self.indentation = ''
        
    def generate_string(self,*strings):
        # incase a tuple, convert to list
        string = ''
        string_list = []
        for entry in strings:
            t = type(entry)
            if t == str:
                s = entry
            elif t == tuple:
                s = ''
                if len(entry) > 0:
                    s = str(entry[0])
            else:
                print_err_msg('Unknown type in Output.generate_string()')
            string_list.append(s)
            string = ''.join(string_list)
        return string
        
    def print_and_return(self,*strings):
        string = self.generate_string(strings)
        if self.fid:
            self.fid.write(self.indentation+string+'\n')
        sys.stdout.write(self.indentation+string+'\n')

    def print_no_return(self,*strings):
        string = self.generate_string(strings)
        if self.fid:
            self.fid.write(self.indentation+string)
        sys.stdout.write(self.indentation+string)

    def destroy(self):
        if self.fid:
            self.fid.close()
