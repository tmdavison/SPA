import numpy as np

def latlon_to_cartesian(latlon):
    # Extract lat, lon from v array
    lat = latlon[0]
    lon = latlon[1]
    
    # Convert lat to colatitude
    colat = 90 - lat
     
    # Convert to radians
    theta = np.radians(colat)
    phi = np.radians(lon)

    # Convert to cartesian coords
    vx = np.sin(theta) * np.cos(phi)
    vy = np.sin(theta) * np.sin(phi)
    vz = np.cos(theta)
    
    return vx, vy, vz

def vector_mag(v):
    v = np.array(v)
    return np.sqrt(v.dot(v))

def get_rotation_axis_and_angle(p1, p2):
    '''
    Given two points on a sphere, find the axis of rotation unit vector and the angle between them
    
    Inputs
    ======
    p1: array-like
        (lat, lon) of starting point
    p2: array-like
        (lat, lon) of final point
        
    Returns
    =======
    krot: array-like
        rotation axis unit vector (cartesian)
    arot: float
        angle of rotation, degrees
    '''
    
    # First convert positions to cartesian vectors

    p1c = latlon_to_cartesian(p1)
    p2c = latlon_to_cartesian(p2)
    
    # Cross product of two vectors gives perpendicular vector (i.e. rotation axis)
    krot = np.cross(p1c, p2c)
    kmag = vector_mag(krot)
    
    # convert to unit vector
    krot /= kmag
    
    # Angle of rotation
    arot = np.degrees(np.arccos(np.dot(p1c, p2c) / (vector_mag(p1c) * vector_mag(p2c))))
    
    return krot, arot

def rotate_coords(v, k, rotangle):
    '''
    Rotates coordinates on a sphere
    
    Inputs
    ------
    v : array-like
        lat, lon of every point to rotate
    k : array-like
        axis of rotation (unit vector)
    rotangle : float
        angle of rotation, degrees
        
    Returns
    -------
    vrotated : array-like
        lat, lon of all rotated points
    '''
        
    # Convert to radians
    rotangle = np.radians(rotangle)

    # Convert to cartesian
    vx, vy, vz = latlon_to_cartesian(v)
    
    # store compondents of each point as an array 
    vcart = np.vstack([vx, vy, vz])

    # extract components of k vector
    kx = k[0]
    ky = k[1]
    kz = k[2]
    
    # ensure k is unit-vector
    kmag = vector_mag(k)
    # kmag should be 1, but if not, we can scale the components here
    kx /= kmag
    ky /= kmag
    kz /= kmag
    
    # Construct K matrix    
    K = np.matrix([[0, -kz, ky],[kz, 0, -kx],[-ky, kx, 0]])
    
    # Find K-squared
    K2 = K * K
    
    # Identity matrix
    I = np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    
    # Compute R matrix
    R = I + (np.sin(rotangle) * K) + ((1 - np.cos(rotangle)) * K2)
    
    # Rotate cartesian coords
    vrot = R * vcart
    
    # Convert vrot to spherical (colatitude, longitude), and convert matrices to arrays (.A1)
    thetarot = np.degrees(np.arccos(vrot[2])).A1
    lonrot = np.degrees(np.arctan2(vrot[1], vrot[0])).A1
    
    # Convert to latitude
    latrot = 90 - thetarot
    
    # Stack lat and lon
    return np.column_stack([latrot, lonrot])

def rotate_coords2(v, p1, p2):
    '''
    Wrapper for rotate coords which takes starting and ending positions to
    calculate the rotation axis and angle and then rotate all the lat,lon pairs
    '''
        
    k, a = get_rotation_axis_and_angle(p1, p2)
    
    return rotate_coords(v, k, a)
