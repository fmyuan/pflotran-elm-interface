import sys
import shutil
import os
import re

swap = re.compile(r'(?P<pre>\s+[0-9]+\.[0-9]*)(?P<post>-[0-9]{3})')

if len(sys.argv) < 2:
    print('A filename must be entered on the command line. E.g.\n\n'
          '    python %s <filename>\n'%sys.argv[0])
    sys.exit(0)

filename = sys.argv[1]
if not os.path.isfile(filename):
    print('File "%s" does not exist within directory "%s".'
          %(filename,os.getcwd()))
    sys.exit(0)

f2 = open(filename+'.tmp','w')
for line in open(filename,'r'):
    f2.write(swap.sub('\g<pre>E\g<post>',line))
f2.close()

shutil.copy(filename+'.tmp',filename)
os.remove(filename+'.tmp')


