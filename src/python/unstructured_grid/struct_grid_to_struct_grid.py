from h5py import *
import numpy
from get_dimensions import get_dimensions

grid_block_a = \
'''GRID
  TYPE STRUCTURED
  NXYZ {} {} {}
  ORIGIN {} {} {}
  DXYZ
'''

grid_block_b = \
'''  /
END
'''

region_bottom = \
'''
REGION bottom
  FACE BOTTOM
  COORDINATES
    -1.d20 -1.d20 0.
     1.d20  1.d20 0.
  /
END
'''

region_top = \
'''
REGION top
  FACE TOP
  COORDINATES
    -1.d20 -1.d20 {}
     1.d20  1.d20 {}
  /
END
'''

region_west = \
'''
REGION west
  FACE WEST
  COORDINATES
    0. -1.d20 -1.d20
    0.  1.d20  1.d20
  /
END
'''

region_east = \
'''
REGION east
  FACE EAST
  COORDINATES
    {} -1.d20 -1.d20
    {}  1.d20  1.d20
  /
END
'''

region_south = \
'''
REGION south
  FACE SOUTH
  COORDINATES
    -1.d20 0. -1.d20
     1.d20 0.  1.d20
  /
END
'''

region_north = \
'''
REGION north
  FACE NORTH
  COORDINATES
    -1.d20 {} -1.d20
     1.d20 {}  1.d20
  /
END
'''

def get_max_coordinate(nx,x0,dx):
    x = x0
    for i in range(nx):
        if len(dx) == 1:
            x += dx[0]
        else:
            x += dx[i]
    return round(x,10)
 

def write_structured_grid(n,d,origin):
    file = open('structured_grid_card.txt','w')
    file.write(grid_block_a.format(n[0],n[1],n[2],
                                   origin[0],origin[1],origin[2]))
    for i in range(3):
        file.write('   ')
        for entry in d[i]:
            file.write(' {}'.format(entry))
        file.write('\n')
    file.write(grid_block_b)
    file.close()
    file = open('structured_region_cards.txt','w')
    file.write(region_west)
    xmax = get_max_coordinate(n[0],origin[0],d[0])
    file.write(region_east.format(xmax,xmax))
    file.write(region_south)
    ymax = get_max_coordinate(n[1],origin[1],d[1])
    file.write(region_north.format(ymax,ymax))
    file.write(region_bottom)
    zmax = get_max_coordinate(n[2],origin[2],d[2])
    file.write(region_top.format(zmax,zmax))
    file.close()

if __name__ == "__main__":

    n, d, origin = get_dimensions()
    write_structured_grid(n,d,origin)

    print('Done generating structured grid.')
