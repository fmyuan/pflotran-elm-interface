# check matrix transpose.py
import sys
import shutil
import os
import fnmatch
import numpy as np

if len(sys.argv) != 3:
    sys.exit("ERROR: 2 keywords: matrix.dat and matrix_transpose.dat in swap_keyword.py")
    
filename1 = sys.argv[1]
filename2 = sys.argv[2]

print('Comparing "{}" and "{}".'.format(filename1,filename2))

f = open(filename1,'r')
lines = f.readlines()
f.close()
num_rows = len(lines)-2
line = lines[2]
words = line.strip().split()
num_cols = len(words)

array1 = np.zeros((num_rows,num_cols))
array2 = np.zeros((num_rows,num_cols))

for irow in range(num_rows):
    words = lines[irow+2].strip().split()
    for icol in range(num_cols):
        array1[irow][icol] = words[icol]
    
# this matrix is transposed
f = open(filename2,'r')
f.readline()
f.readline()
for icol in range(num_cols):
    words = f.readline().strip().split()
    for irow in range(num_rows):
        array2[irow][icol] = words[irow]
f.close()

max_rel_diff = 0.
for irow in range(num_rows):
    for icol in range(num_cols):
        value1 = array1[irow][icol]
        value2 = array2[irow][icol]
        print(irow+1,icol+1,value1,value2)
        if abs(value1) > 0.:
            max_rel_diff = max(max_rel_diff,abs((value1-value2)/value1))
        elif abs(value2) > 0.:
            print('Both values are not zero.')
            
print('Maximum difference: {}'.format(max_rel_diff))
            
          