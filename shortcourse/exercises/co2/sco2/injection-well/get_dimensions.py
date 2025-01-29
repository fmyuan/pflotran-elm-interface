def get_dimensions():

    # for 543 irregular grid spacing, dx,dy,dz must be specified.
    nx = 20; ny = 1; nz = 25
    dx = [5.15874, 6.70684968, 8.71862112, 11.3342014, 14.7344892, 19.1548207, 24.9012456, 32.3716498, 42.0831569, 54.7080643, 71.120508, 92.4566299, 120.193643, 156.251758, 203.127254, 264.065431, 343.285084, 446.270604, 580.151748, 551.2055014]
    dy = [1e3]
    dz = [24.9936, 18.28800, 9.14400, 11.8872, 11.8872, 9.4488, 7.0104, 7.3152, 6.7056, 6.7056, 7.0104, 6.4008, 8.5344, 6.09600, 3.9624, 2.1336, 7.9248, 3.6576, 3.04800, 1.2192, 4.2672, 3.04800, 3.04800, 3.04800, 15.24]

    # for uniform grid spacing, specify one value in each direction
    #nx = 20; ny = 1; nz = 25
    #dx = [50.]
    #dy = [50.]
    #dz = [10.]

    # note that you can mix and match irregular spacing along the different axes
    n = (nx,ny,nz)
    d = [dx,dy,dz]
    origin = (0.,0.,-1157.024)

    return n, d, origin
