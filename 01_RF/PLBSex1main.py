'''PYTHON WORKFLOW'''

# Import of the required libraries
import numpy as np
import math
from functions import *
from math import sqrt
'''
REMEMBER: in naming variables:
   - P is for Coordinate Vector
   - C is for Covariance Matrix
   - V is for Velocity Vector
'''
'''Input data'''
# ITRF baseline	from COMO to BRUN
P_ComoBrun_GC = np.array([[-1040.168], [-72.970], [1631.398]])
C_ComoBrun_GC = np.array([ [2, 0.5, 0.5],
	                       [0.5, 1, 0.5],
	                       [0.5, 0.5, 2]
	                     ])*10**(-6)
# ITRF baseline between Brunate and Point 1 
P_Brun001_GC = np.array([[-51.130], [76.749], [38.681]])
C_Brun001_GC = np.array([ [1.5, 0.3, 0.3],
	                      [0.3, 1, 0.2],
	                      [0.3, 0.2, 2]
	                    ])*10**(-6)
# Brunate coordinates in LL
P_Brun_LL = np.array([[0], [0], [0]])
# Brunate vertical deflection components
qsi=10.23   #seconds
eta=9.5     #seconds
# Points 1, 2, 3 in Local Level survayed with respect to Brunate
P_1_LL = np.array([[99.9994], [0.0], [1.0174]])
P_2_LL = np.array([[80.0], [85.0], [2.0]])
P_3_LL = np.array([[-25.0],	[10.0],	[-3.0]])
# Covariance matrix equal for the three points in LL
C_LL = np.array([[0.1**2, 0, 0], [0, 0.1**2, 0], [0, 0, 0.15**2]])

'''
STEP 1:  definition of Como Geocentric Cartesian Coordinates for a specific epoch t = 02/02/2019
'''
t_0 = 2010
t = 2019 + (31 + 2)/365
delta_t = t - t_0
# Retrieve from EPN website Como GC coordinates and velocities at t_0
P_Como_GC_2010 = np.array([[4398306.209],[704149.948],[4550154.733]])
V_Como_GC_2010 = np.array([[-0.0146],[0.0180],[0.0112]])
#Como GC at t
P_Como_GC=P_Como_GC_2010+(delta_t*V_Como_GC_2010)

'''STEP 2: computation of Brunate GC'''
P_Brun_GC=P_Como_GC+P_ComoBrun_GC

'''STEP 3: Compute Local Cartesian coord of P1 w.r.t. Brun'''
# Compute P1 in GC (remember: ITRF baseline between Brunate and Point 1 is known)
P_1_GC=P_Brun_GC+P_Brun001_GC
# Transform P_1_GC to Geodetic (use cartToGeod function from functions module)
P1_geod=cartToGeod(P_1_GC)
# Transform P_1 to LC with respect to Brunate (use GCtoLC)
x,y,z=P1_geod[0],P1_geod[1],P1_geod[2]
P_Brun_GG = cartToGeod(P_Brun_GC)
phiBrun, lambdaBrun = degToRad(P_Brun_GG[0][0]), degToRad(P_Brun_GG[1][0])
Brun1_LC = GCtoLC(P_Brun001_GC, phiBrun, lambdaBrun)
P_Brun1_LC = Brun1_LC[0]

'''STEP 4: Convert LC of 0001 to pseudo Local Level w.r.t. Brun'''
# Define Rx and Ry matrix and compute the temporary value required to get alpha
ksi=np.array([0,0,10.23])
eta=np.array([0,0,9.5])
degksi=sexToDeg(ksi)
degeta=sexToDeg(eta)
radksi=degToRad(degksi)
radeta=degToRad(degeta)
#Rx e Ry
Rx=np.array([[1,0,0],[0,1,-radksi],[0,radksi,1]])
Ry=np.array([[1,0,-radeta],[0,1,0],[radeta,0,1]])
P_Brun1_temp=np.dot((np.dot(Ry,Rx)),P_Brun1_LC)

'''STEP 5: Compute alpha rotation between LC and LL in Brun'''
# Compute alpha and Rz
radalfa=math.atan(P_Brun1_temp[1]/P_Brun1_temp[0])
degalfa=radToDeg(radalfa)
sexalfa=degToSex(degalfa)
#Rz
cosalfa=math.cos(radalfa)
sinalfa=math.sin(radalfa)
Rz=np.array([[cosalfa,sinalfa,0],[-sinalfa,cosalfa,0],[0,0,1]])
#LL of 0001 wrt Brunate
P_Brun1_LL=np.dot(np.dot(np.dot(Rz,Ry),Rx),P_Brun1_LC)

'''STEP 6: Compute LC of P2 and P3'''
R=np.dot(np.dot(Rz,Ry),Rx)
Rt=np.transpose(R)
P_Brun2_LC=np.dot(Rt,P_2_LL)
P_Brun3_LC=np.dot(Rt,P_3_LL)

'''STEP 7: Compute GC of P2 and P3'''
sinphi=math.sin(phiBrun)
sinlambda=math.sin(lambdaBrun)
cosphi=math.cos(phiBrun)
coslambda=math.cos(lambdaBrun)
R0=np.array([[-sinlambda,coslambda,0],[(-sinphi*coslambda),(-sinphi*sinlambda),cosphi],[(cosphi*coslambda),(cosphi*sinlambda),sinphi]])
R0t=np.transpose(R0)
baseline2=np.dot(R0t,P_Brun2_LC)
baseline3=np.dot(R0t,P_Brun3_LC)
P_2_GC=baseline2+P_Brun_GC
P_3_GC=baseline3+P_Brun_GC

'''STEP 8: Conversion through EPN website ITRS GC coord of the stations (Brun, P1, P2, P3) to ETRF'''
P_Como_ETRF = [4398306.50780, 704149.56120, 4550154.50290]      
P_Brun_ETRF = [4397266.33990, 704076.59130, 4551785.90090]      
P_1_ETRF = [4397215.20990, 704153.34030, 4551824.58190]      
P_2_ETRF = [4397183.08940, 704084.32790, 4551867.37790]      
P_3_ETRF = [4397272.15520, 704050.77260, 4551780.10640]
#conversion in ETRF geodetic
P_Como_ETRF_GG = newcartToGeod(P_Como_ETRF)
P_Brun_ETRF_GG = newcartToGeod(P_Brun_ETRF)
P_1_ETRF_GG = newcartToGeod(P_1_ETRF)
P_2_ETRF_GG = newcartToGeod(P_2_ETRF)
P_3_ETRF_GG = newcartToGeod(P_3_ETRF)

''' STEP 9: covariance matrixes'''
#ITRF standard deviations and velocity for como
C_Como_std_t0 = np.array([(0.001**2, 0, 0), (0, 0.001**2, 0), (0, 0, 0.001**2)])
C_Como_V_t0 = np.array([(0.0001**2, 0, 0), (0, 0.0001**2, 0), (0, 0, 0.0001**2)])
#covariances in GC
C_Como_t = np.array(C_Como_std_t0 + np.dot(delta_t**2, C_Como_V_t0))
C_Brun_t = C_Como_t + C_ComoBrun_GC
C_001_t = C_Brun_t + C_Brun001_GC # Geocentric covariance
#covariances in LC
C1_LC = np.dot(np.dot(R0, C_001_t), R0.transpose())
C_LC_Brunate = np.dot(np.dot(R0,C_Brun_t),R0.transpose())
C_LC=np.dot(np.dot(R0.transpose(), C_LL), R0)
C_LC_2_3= C_LC + C_LC_Brunate
C_LC_2_3 = np.dot(np.dot(R0,C_LC_2_3),R0.transpose())

''' STEP 10: standard deviations''' 
Brunate_std = sqrt(C_LC_Brunate[0][0]), sqrt(C_LC_Brunate[1][1]), sqrt(C_LC_Brunate[2][2])
C1_LC_std = sqrt(C1_LC[0][0]), sqrt(C1_LC[1][1]), sqrt(C1_LC[2][2])
C_LC_2_3_std = sqrt(C_LC_2_3[0][0]), sqrt(C_LC_2_3[1][1]), sqrt(C_LC_2_3[2][2])

#printing results
print('COMO in ITRF at present epoch: \n', P_Como_GC)
print('BRUNATE in ITRF at present epoch: \n', P_Brun_GC)
print('P1 in Local Cartesian \n', P_Brun1_LC)
print('P2 in Local Cartesian: \n', P_Brun2_LC)
print('P3 in Local Cartesian: \n', P_Brun3_LC)
print('Brunate in Geodetic:\n',P_Brun_GG)
print('phi e lambda in radiant:',phiBrun,lambdaBrun)
print('Alpha in sexadecimal: \n', sexalfa)
print('P1 in ITRF Geocentric cartesian: \n', P_1_GC)
print('P2 in ITRF Geocentric Cartesian: \n', P_2_GC)
print('P3 in ITRF Geocentric Cartesian: \n', P_3_GC)
print('P1 in Geodetic:\n',P1_geod)
print('P1 in Local Level:',P_Brun1_LL)
print('Covariance matrix for P1 in LC: \n', C1_LC)
print('Covariance matrix for P2 e P3 in LC: \n', C_LC_2_3)
print('Covariance matrix for Brunate in LC : \n', C_LC_Brunate)
print('Como in ETRF: \n', P_Brun_ETRF)
print('Brunate in ETRF: \n', P_Brun_ETRF)
print('P1 in ETRF: \n', P_1_ETRF)
print('P2 in ETRF: \n', P_2_ETRF)
print('P3 in ETRF: \n', P_3_ETRF)
print('Como in ETRF Geodetic: \n', P_Como_ETRF_GG)
print('Brun in ETRF Geodetic : \n', P_Brun_ETRF_GG)
print('P1 in ETRF Geodetic: \n', P_1_ETRF_GG)
print('P2 in ETRF Geodetic: \n', P_2_ETRF_GG)
print('P3 in ETRF Geodetic: \n', P_3_ETRF_GG)
print('std of Brunate in [m]:\n', Brunate_std)
print('std of P1 in [m]:\n', C1_LC_std)
print('std of P2 and P3 in [m]:\n', C_LC_2_3_std)



