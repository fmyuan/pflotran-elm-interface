# diff_matrices.py
import sys
import numpy as np

if len(sys.argv) < 3:
    sys.exit("ERROR: Must include at least two PETSc vector files")
    
filename1 = sys.argv[1]
filename2 = sys.argv[2]

verbose = False
if len(sys.argv) == 4:
    word = sys.argv[3]
    verbose = (word.lower() == 'verbose')
    if not transpose:
        sys.exit("ERROR: Unrecognized flag. Did you mispell 'verbose'?")

print('Comparing "{}" and "{}".'.format(filename1,filename2))

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

f = open(filename1,'r')
lines = []
for line in f:
    if is_number(line):
        lines.append(line)
f.close()
num_rows = len(lines)

array1 = np.zeros(num_rows)
array2 = np.zeros(num_rows)

i = 0
for line in lines:
    array1[i] = float(line)
    i += 1

f = open(filename2,'r')
i = 0
for line in f:
    if is_number(line):
        array2[i] = float(line)
        i += 1
f.close()

if not i == num_rows:
    sys.exit("ERROR: Number of entries in vectors does not match ({} vs {})".format(i,num_rows))

max_abs_val = max(np.absolute(array1).max(),np.absolute(array2).max())
min_abs_val = min(np.absolute(array1).min(),np.absolute(array2).min())
max_rel_diff = 0.
tpl = (0.,0.,0.,0.,0)
num_max = 10
max_tpl = [tpl]*num_max
for irow in range(num_rows):
    value1 = array1[irow]
    value2 = array2[irow]
    rel_diff = 0.
    diff = 0.
    diff = value1-value2
    if abs(value1) > 0.:
        rel_diff = diff/value1
    elif abs(value2) > 0.:
        rel_diff = diff/value2
    if verbose:
        print('{} {:.8e} {:.8e} {:.5f} {:.3e}'.format(irow+1,value1,value2,rel_diff,diff/max_abs_val))
    for i in range(num_max):
        if abs(rel_diff) > abs(max_tpl[i][0]):
            max_tpl[i] = (rel_diff,value1,value2,irow+1)
            break
        
for i in range(num_max):
    tpl = max_tpl[i]
    print('{} Maximum relative difference: {:.3e} (values = {:.8e}, {:.8e}) at {}'.format(i,tpl[0],tpl[1],tpl[2],tpl[3]))
print('Maximum absolute value: {:.8e}'.format(max_abs_val))
print('Minimum absolute value: {:.8e}'.format(min_abs_val))
