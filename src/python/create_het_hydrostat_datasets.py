'''
  Checkout glenn/hydrostatic-printout to print out the hydrostatic_*.txt files.
  Then, run through script to generate dataset.
'''
import sys
from h5py import *
import numpy

def read_hydrostatic_file(filename):
  f = open(filename,'r')
  array = []
  for line in f:
    w = line.split()
    l = []
    l.append(int(w[2]))
    l.append(int(w[1]))
    l.append(float(w[3]))
    array.append(l)
  return array

def create_map(direction,array):
  count = 0
  for i in range(len(array)):
    if array[i][0] == direction:
      count += 1
  map_array = numpy.zeros((count,2),'=i4')
  count = 0
  for i in range(len(array)):
    if array[i][0] == direction:
      map_array[count][0] = count+1
      map_array[count][1] = array[i][1]
      count += 1
  return map_array 
  
def create_dataset(direction,arrays):
  count = 0
  for icell in range(len(arrays[0])):
    if arrays[0][icell][0] == direction:
      count += 1
  data_array = numpy.zeros((count,len(arrays)),'=f8')
  for iarray in range(len(arrays)):
    count = 0
    for icell in range(len(arrays[iarray])):
      if arrays[iarray][icell][0] == direction:
        data_array[count][iarray] = arrays[iarray][icell][2]
        count += 1
  return data_array 

def create_region(direction,array):
  count = 0
  for i in range(len(array)):
    if array[i][0] == direction:
      count += 1
  cell_id_array = numpy.zeros(count,'=i4')
  face_id_array = numpy.zeros(count,'=i4')
  count = 0
  for i in range(len(array)):
    if array[i][0] == direction:
      cell_id_array[count] = array[i][1]
      face_id_array[count] = array[i][0]
      count += 1
  return cell_id_array, face_id_array

filenames = []
filenames.append('hydrostatic_0.000000000000000E+00.txt')
filenames.append('hydrostatic_1.000000000000000E+01.txt')
filenames.append('hydrostatic_5.000000000000000E+01.txt')
filenames.append('hydrostatic_1.000000000000000E+02.txt')
'''
filenames.append('hydrostatic_0.000000000000000E+00.txt')
filenames.append('hydrostatic_1.000000000000000E-06.txt')
filenames.append('hydrostatic_3.000000000000000E-06.txt')
filenames.append('hydrostatic_7.000000000000000E-06.txt')
filenames.append('hydrostatic_1.500000000000000E-05.txt')
filenames.append('hydrostatic_3.100000000000000E-05.txt')
filenames.append('hydrostatic_6.299999999999999E-05.txt')
filenames.append('hydrostatic_1.270000000000000E-04.txt')
filenames.append('hydrostatic_2.550000000000000E-04.txt')
filenames.append('hydrostatic_5.109999999999998E-04.txt')
filenames.append('hydrostatic_1.023000000000000E-03.txt')
filenames.append('hydrostatic_2.047000000000000E-03.txt')
filenames.append('hydrostatic_4.095000000000000E-03.txt')
filenames.append('hydrostatic_8.190999999999999E-03.txt')
filenames.append('hydrostatic_1.638300000000000E-02.txt')
filenames.append('hydrostatic_3.276700000000000E-02.txt')
filenames.append('hydrostatic_6.553500000000000E-02.txt')
filenames.append('hydrostatic_1.310710000000000E-01.txt')
filenames.append('hydrostatic_2.621430000000000E-01.txt')
filenames.append('hydrostatic_5.242870000000001E-01.txt')
filenames.append('hydrostatic_1.048575000000000E+00.txt')
filenames.append('hydrostatic_2.097151000000000E+00.txt')
filenames.append('hydrostatic_4.194303000000000E+00.txt')
filenames.append('hydrostatic_8.388606999999999E+00.txt')
filenames.append('hydrostatic_1.677721499999999E+01.txt')
filenames.append('hydrostatic_2.677721500000000E+01.txt')
filenames.append('hydrostatic_3.677721500000000E+01.txt')
filenames.append('hydrostatic_4.677721499999999E+01.txt')
filenames.append('hydrostatic_5.677721500000001E+01.txt')
filenames.append('hydrostatic_6.677721500000000E+01.txt')
filenames.append('hydrostatic_7.677721500000000E+01.txt')
filenames.append('hydrostatic_8.677721500000000E+01.txt')
filenames.append('hydrostatic_9.677721500000000E+01.txt')
filenames.append('hydrostatic_1.000000000000000E+02.txt')
'''

arrays = []
for filename in filenames:
  arrays.append(read_hydrostatic_file(filename))

times_array = numpy.zeros(len(arrays),'=f8')
for i in range(len(times_array)):
  w = filenames[i].split('_')
  w = w[1].split('.t')
  times_array[i] = float(w[0])

directions = []
directions.append(2)
directions.append(3)
directions.append(4)
directions.append(6)

filename = 'dataset.h5'
h5file = File(filename,mode='w')

for i in range(len(directions)):
  direction = directions[i]
  grp_name = 'Map_%d' % direction
  m = create_map(direction,arrays[0])
  h5file[grp_name + '/Data'] = m
  grp_name = 'Pressure_%d' % direction
  d = create_dataset(direction,arrays)
  h5file[grp_name + '/Data'] = d
  h5file[grp_name + '/Times'] = times_array
  h5grp = h5file[grp_name]
  h5grp.attrs['Time Units'] = numpy.string_('d')

h5file.close()

filename = 'regions.h5'
h5file = File(filename,mode='w')

for i in range(len(directions)):
  direction = directions[i]
  grp_name = 'east_%d' % direction
  cell_ids, face_ids = create_region(direction,arrays[0])
  h5file['Regions/' + grp_name + '/Cell Ids'] = cell_ids
  h5file['Regions/' + grp_name + '/Face Ids'] = face_ids

h5file.close()
