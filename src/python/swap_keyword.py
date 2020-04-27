# swap_keyword.py
import sys
import shutil
import os
import fnmatch

if len(sys.argv) != 3:
    sys.exit("ERROR: 2 keywords: old new must be supplied in swap_keyword.py")

keyword_old = sys.argv[1]
keyword_new = sys.argv[2]

def swap_file(filename):
    f = open(filename,'r')
    f2 = open(filename+'.tmp','w')
    for line in f:
        line2 = line.replace(keyword_old,keyword_new)
        f2.write(line2)
    f.close()
    f2.close()
    # using shutil.move adds ^M to end of lines.
    os.remove(filename)
    shutil.copy(filename+'.tmp',filename)
    os.remove(filename+'.tmp')

suffix = '*.in'
for root, dirnames, filenames in os.walk('.'):
    for filename in fnmatch.filter(filenames,suffix):
        filename = os.path.join(root,filename)
        print(filename)
        swap_file(filename)

print('done')
