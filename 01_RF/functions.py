import math
import numpy as np

def radToDeg(rad):
    deg = rad*180/math.pi
    return deg

def degToRad(deg):
    rad = deg*math.pi/180
    return rad

def sexToDeg(sex):
	# remember sex = [deg, min, sec]
	deg = sex[0] + sex[1]/60 + sex[2]/3600
	return deg

def degToSex(d):
	deg = np.fix(d)
	minu =(d-deg)*60
	sec =(minu - np.fix(minu))*60
	return [deg, minu, sec]

def cartToGeod(P_cart):
    '''
    Input:
       - P_cart = numpy array containing cartesian coordinates of the point
    '''
    x, y, z = P_cart[0][0], P_cart[1][0], P_cart[2][0]
    a = 6378137
    e = 0.0818191908426215
    e2 = e*e
    b = a*(np.sqrt(1 - e2))
    eb2 = (a*a - b*b)/(b*b)
    # radius computation 
    r = np.sqrt(x*x + y*y)
    # longitude
    lon = np.arctan2(y, x)
    #latitute
    psi = np.arctan2(z, (r*np.sqrt(1-e2)))
    lat = np.arctan2((z+eb2*b*np.power(np.sin(psi), 3)), (r - e2*a*np.power(np.cos(psi), 3)))
    N = a/np.sqrt(1 - e2*np.power(np.sin(lat), 2))
    h = r/(np.cos(lat)) - N
    lon = radToDeg(lon)
    lat = radToDeg(lat)    

    return np.array([[lat], [lon], [h]])

def newcartToGeod(P_cart):
    '''
    Input:
       - P_cart = numpy array containing cartesian coordinates of the point
    '''
    x, y, z = P_cart[0], P_cart[1], P_cart[2]
	#a = 6378137
    a = 6378137
    e = 0.0818191908426215
    e2 = e*e
    b = a*(np.sqrt(1 - e2))
    eb2 = (a*a - b*b)/(b*b)
    # radius computation 
    r = np.sqrt(x*x + y*y)
    # longitude
    lon = np.arctan2(y, x)
    #latitute
    psi = np.arctan2(z, (r*np.sqrt(1-e2)))
    lat = np.arctan2((z+eb2*b*np.power(np.sin(psi), 3)), (r - e2*a*np.power(np.cos(psi), 3)))
    N = a/np.sqrt(1 - e2*np.power(np.sin(lat), 2))
    h = r/(np.cos(lat)) - N
    lon = radToDeg(lon)
    lat = radToDeg(lat)
    return np.array([[lat], [lon], [h]])

def geodToCart(P_geod):
    lat, lon, h = P_geod[0][0], P_geod[1][0], P_geod[2][0]
    a = 6378137
    e = 0.0818191908426215
    e2 = e*e
    lat = degToRad(lat)
    lon = degToRad(lon)
    N = a/np.sqrt(1 - e2*np.power(np.sin(lat), 2))
    x = (N + h)*np.cos(lat)*np.cos(lon)
    y = (N + h)*np.cos(lat)*np.sin(lon)
    z = (N*(1-e2) + h)*np.sin(lat)
    
    return np.array([[x], [y], [z]])

def GCtoLC(P_GC, phi, lambdaValue):
    '''
    Input:
       - P_GC = baseline between the point and the center of LC
       - phi = latitude of the center of LC in radiants
       - lambda = longitude of the center of LC in radiants
    '''
    # Definition of the rotation matrix R0

    R0 = np.array([ [-np.sin(lambdaValue), np.cos(lambdaValue), 0],
                    [-np.sin(phi)*np.cos(lambdaValue), -np.sin(phi)*np.sin(lambdaValue), np.cos(phi)],
                    [np.cos(phi)*np.cos(lambdaValue), np.cos(phi)*np.sin(lambdaValue), np.sin(phi)]
                  ])
    P_LC = np.dot(R0, P_GC)
    #C_LC = np.dot(R0, C_GC, R0)

    return [P_LC]