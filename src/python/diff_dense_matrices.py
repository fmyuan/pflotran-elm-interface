# diff_matrices.py
import sys
import numpy as np

if len(sys.argv) < 3:
    sys.exit("ERROR: Must include at least two matrix files and a potential 'transpose' flag in diff_matrix.py")
    
filename1 = sys.argv[1]
filename2 = sys.argv[2]
transpose = False
if len(sys.argv) == 4:
    word = sys.argv[3]
    transpose = (word.lower() == 'transpose')
    if not transpose:
        sys.exit("ERROR: Unrecognized flag. Did you mispell 'tranpose'?")
    

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

num_rows2 = num_rows
num_cols2 = num_cols
if transpose:
    print('transposing second matrix')
    num_rows2 = num_cols
    num_cols2 = num_rows
    
f = open(filename2,'r')
f.readline()
f.readline()
row_count = 0
for irow in range(num_rows2):
    words = f.readline().strip().split()
    if not len(words) == num_cols2:
        i = len(words)
        if transpose:
            sys.exit("ERROR: Number of columns in second matrix does not match the number of rows in the first ({} vs {})".format(i,num_cols2))
        else:
            sys.exit("ERROR: Number of columns in second matrix does not match the number of columns in the first ({} vs {})".format(i,num_cols2))
    for icol in range(num_cols2):
        if transpose:
            array2[icol][irow] = words[icol]
        else:
            array2[irow][icol] = words[icol]
    row_count += 1
f.close()

if not row_count == num_rows2:
    if transpose:
        sys.exit("ERROR: Number of rows in second matrix does not match the number of columns in the first ({} vs {})".format(row_count,num_rows2))
    else:
        sys.exit("ERROR: Number of rows in second matrix does not match the number of rows in the first ({} vs {})".format(row_count,num_rows2))

max_abs_val = max(np.absolute(array1).max(),np.absolute(array2).max())
min_abs_val = min(np.absolute(array1).min(),np.absolute(array2).min())
max_rel_diff = 0.
tpl = (0.,0.,0.,0.,0,0)
num_max = 10
max_tpl = [tpl]*num_max
for irow in range(num_rows):
    for icol in range(num_cols):
        value1 = array1[irow][icol]
        value2 = array2[irow][icol]
        rel_diff = 0.
        diff = 0.
        if abs(value1) > 0.:
            diff = value1-value2
            rel_diff = diff/value1
        elif abs(value2) > 0.:
            print('Both values are not zero.')
        print('{} {} {:.8e} {:.8e} {:.5f} {:.3e}'.format(irow+1,icol+1,value1,value2,rel_diff,diff/max_abs_val))
        for i in range(num_max):
            if abs(rel_diff) > abs(max_tpl[i][0]):
                 max_tpl[i] = (rel_diff,value1,value2,irow+1,icol+1)
                 break
            
for i in range(num_max):
    tpl = max_tpl[i]
    print('{} Maximum relative difference: {:.3e} (values = {:.8e}, {:.8e}) at {} {}'.format(i,tpl[0],tpl[1],tpl[2],tpl[3],tpl[4]))
print('Maximum absolute value: {:.8e}'.format(max_abs_val))
print('Minimum absolute value: {:.8e}'.format(min_abs_val))
