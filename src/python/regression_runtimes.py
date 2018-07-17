# list regression test run times in descending order
import sys

times = []
for line in open(sys.argv[1],'r'):
  if line.find('run time :') > -1:
    times.append(line.strip())

flag = True
while(flag):
  flag = False
  for i in range(len(times)-1):
    if float(times[i].split(':')[2].split()[0]) < \
       float(times[i+1].split(':')[2].split()[0]):
      times[i], times[i+1] = times[i+1],times[i]
      flag = True

for i in range(len(times)):
  print(times[i])

print('done')
