"""
Description: Error messaging for python scripts
Author: Glenn Hammond
"""

import sys

def print_err_msg(*strings):
    list = []
    for string in strings:
        list.append(string)
    string = ''.join(list)
    print(string)
    sys.exit(1)