"""
Description: Defines a progress bar class for printing progress of an iterative
             process of known total number of iterations to the screen.
Author: Glenn Hammond
"""
import sys

class ProgressBar():
    
    def __init__(self,maximum_iteration,print_fraction=0.1):
        self.maximum_iteration = maximum_iteration
        self.print_fraction = print_fraction
        self.print_increment = int(self.print_fraction*self.maximum_iteration)
        self.iteration = 0
        self.print_flag = True
        if self.maximum_iteration < 100:
            self.print_flag = False
        if self.print_flag:
            # in the future, use print("",end='') avialable with python3
            sys.stdout.write('Progress: 0%')
            sys.stdout.flush()
        
    def increment(self):
        if not self.print_flag:
             return 
        self.iteration += 1
        if self.iteration == self.maximum_iteration or \
           self.iteration % self.print_increment == 0:
            self.print_progress()
            sys.stdout.flush()
            
    def print_progress(self):
        if self.iteration == self.maximum_iteration:
            sys.stdout.write('-100%\n')
        else:
            percent = int(int((1.e-4+
                               float(self.iteration) / self.maximum_iteration) /
                              self.print_fraction)*self.print_fraction*100.)
            sys.stdout.write('-{}%'.format(percent))
            if percent == 100:
                self.print_flag = False
                sys.stdout.write('\n')
        
        
        
    
