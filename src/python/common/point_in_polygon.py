"""
Description: Defines a function that calculates whether a point lies within
             a prescribed polygon with at least two (rectangle) or more 
             (polygon) points.
Author: Glenn Hammond
"""
import sys

def PointInPolygon(point,polygon,orientation='XY',debug=False):
    """
    point       : [x,y] or [x,y,z]
    polygon     : [[x0,y0],[x1,y1],...,[xn,yn]]
    orientation : string - 'XY','XZ','YZ'
    """
    num_points_in_polygon = len(polygon)
    x = -1.e20
    y = -1.e20
    if orientation == 'XY':
         x = point[0]
         y = point[1]
    elif orientation == 'XZ':
         x = point[0]
         y = point[2]
    elif orientation == 'YZ':
         x = point[1]
         y = point[2]
    else:
        sys.exit('Unrecognized orientation in PointInPolygon(): {}'.
                 format(orientation))
    if num_points_in_polygon == 2:
        x0 = min(polygon[0][0],polygon[1][0])
        x1 = max(polygon[0][0],polygon[1][0])
        y0 = min(polygon[0][1],polygon[1][1])
        y1 = max(polygon[0][1],polygon[1][1])
        in_polygon = (x >= x0 and x <= x1 and y >= y0 and y <= y1)
    else:
        in_polygon = False
        i = num_points_in_polygon - 1
        for j in range(num_points_in_polygon):
            if (polygon[j][1] < y and polygon[i][1] >= y) or \
               (polygon[i][1] < y and polygon[j][1] >= y):
                if (polygon[j][0] + 
                    (y-polygon[j][1])/(polygon[i][1]-polygon[j][1])*
                    (polygon[i][0]-polygon[j][0])) < x:
                    in_polygon = not in_polygon
            i = j
    if debug:
        point_string = 'Point <{},{}>'.format(x,y)
        if in_polygon:
            print('{} is in polygon.'.format(point_string))
        else:
            print('{} is NOT in polygon.'.format(point_string))
    return in_polygon        
        
def TestPointInPolygon():
    orientation = 'YZ'
    debug = True

    polygon = []
    polygon.append([0.,0.])
    polygon.append([1.,0.])
    polygon.append([1.,1.])
    polygon.append([0.,1.])

    point0 = [53235,0.5,0.5]
    point1 = [1,1.000001,0.5]
    point2 = [-1e50,-.000001,0.5]
    point3 = [0.,0.5,1.000001]
    point4 = [1e50,0.5,-.0000001]

    PointInPolygon(point0,polygon,orientation,debug)
    PointInPolygon(point1,polygon,orientation,debug)
    PointInPolygon(point2,polygon,orientation,debug)
    PointInPolygon(point3,polygon,orientation,debug)
    PointInPolygon(point4,polygon,orientation,debug)

    print('')
    polygon = []
    polygon.append([0.,0.])
    polygon.append([0.,1.])
    polygon.append([1.,1.])
    polygon.append([1.,0.])

    PointInPolygon(point0,polygon,orientation,debug)
    PointInPolygon(point1,polygon,orientation,debug)
    PointInPolygon(point2,polygon,orientation,debug)
    PointInPolygon(point3,polygon,orientation,debug)
    PointInPolygon(point4,polygon,orientation,debug)

    print('')
    polygon = []
    polygon.append([1.,1.])
    polygon.append([0.,0.])

    PointInPolygon(point0,polygon,orientation,debug)
    PointInPolygon(point1,polygon,orientation,debug)
    PointInPolygon(point2,polygon,orientation,debug)
    PointInPolygon(point3,polygon,orientation,debug)
    PointInPolygon(point4,polygon,orientation,debug)    

if __name__ == "__main__":
    try:
        TestPointInPolygon()
    except Exception as error:
        print(str(error))
#        if cmdl_options.backtrace:
#            traceback.print_exc()
#        traceback.print_exc()
        print("failure")
        sys.exit(1)

