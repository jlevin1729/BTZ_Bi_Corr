#!/usr/bin/env python3
import numpy as np
import math as m
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as pat
import scipy.optimize as opt
import random as rand
import pickle
#sys.path.append('/home/jole1266/BTZ_geod.py')

#Define a whole bunch of functions. This should only have to be run once.

#Coordinate transforms between Cartesian and polar:

def PolarToCartesian(vec):
    #Input a 2-tuple in polar coordinates (r,theta), get 2-tuple (x,y)
    return(vec[0]*m.cos(vec[1]),vec[0]*m.sin(vec[1]))

def CartesianToPolar(vec):
    #Input a 2-tuple in Cartesian coordinates (x,y), get 2-tuple (r,theta)
    return(m.sqrt(vec[0]**2+vec[1]**2),np.arctan2(vec[1],vec[0]))

#Coordinate transformations between BTZ disk, BTZ geometry, and Poincare disk: 

def BTZDiskToBTZ(r):
    #Input radial coord between 0 and 1
    return(2*r/(1 - r**2))

def BTZToBTZDisk(r):
    return((m.sqrt(1 + r**2)-1)/r)

def BTZToPoincareDisk(vec, rH):
    #Input a tuple (r,phi) in BTZ coordinates along with BTZ horizon radius rH
	
    B = (vec[0]*m.cosh(rH*vec[1]) - rH)/(vec[0]*m.cosh(rH*vec[1]) + rH)

    if B < 0:

        B = 0

    rPoin = m.sqrt(B)

    if abs(vec[0] - rH) < 0.0000000001:

        if vec[1] >= 0.0:

            phiPoin = m.pi/2

        else:

            phiPoin = -m.pi/2

    else:

        if vec[1] >= 0.0:

            A = ((vec[0])**2 - (rH)**2)/((vec[0])**2 * (m.cosh(rH*vec[1]))**2 - (rH)**2)

            if A < 0:
	
                A = 0
           
            phiPoin = m.acos(m.sqrt(A))
	
        else:
		
            A = ((vec[0])**2 - (rH)**2)/((vec[0])**2 * (m.cosh(rH*vec[1]))**2 - (rH)**2)

            if A < 0:
	
                A = 0

            phiPoin = -m.acos(m.sqrt(A))

    return(rPoin, phiPoin)

def PoincareDiskToBTZ(vec, rH):
    rBTZ = rH * m.sqrt(1 + 4*(vec[0])**2 * (m.cos(vec[1]))**2/(1-vec[0]**2)**2)
    if vec[1] >= 0.0:
        phiBTZ = np.arccosh(m.sqrt((1+vec[0]**2)**2/((1-vec[0]**2)**2 + 4 * vec[0]**2 *(m.cos(vec[1]))**2)))/rH
    else:
        phiBTZ = -np.arccosh(m.sqrt((1+vec[0]**2)**2/((1-vec[0]**2)**2 + 4 * vec[0]**2 *(m.cos(vec[1]))**2)))/rH
    return(rBTZ, phiBTZ)

def BTZDiskToPoincareDisk(vec, rH):
    #Input a tuple (r, phi) on the BTZ disk, transform the radial cooordinate to BTZ and then go to Poincare Disk
    BTZvec = (BTZDiskToBTZ(vec[0]),vec[1])
    return(BTZToPoincareDisk(BTZvec, rH))

def PoincareDiskToBTZDisk(vec, rH):
    BTZvec = PoincareDiskToBTZ(vec, rH)
    return(BTZToBTZDisk(BTZvec[0]), BTZvec[1])

#The data of the circle making the geodesic between two points on the Poincare disk, and the geodesic length:

def PoincareGeodesicCircle(vec1, vec2):
    delta_phi = np.arctan2(vec2[1],vec2[0]) - np.arctan2(vec1[1],vec1[0])
    r1 = m.sqrt(vec1[0]**2 + vec1[1]**2)
    r2 = m.sqrt(vec2[0]**2 + vec2[1]**2)
    c_x = -(vec1[1]*(1 + r2**2) - vec2[1]*(1 + r1**2))/(2*r1*r2*m.sin(delta_phi))
    c_y = (vec1[0]*(1 + r2**2) - vec2[0]*(1 + r1**2))/(2*r1*r2*m.sin(delta_phi))
    R = m.sqrt((c_x - vec1[0])**2 + (c_y - vec1[1])**2)
    return((c_x,c_y),R)

def AdS_distance(vec): #two AdS cartesian points
    vec1 = vec[0]
    vec2 = vec[1]
    Deltax = vec2[0]-vec1[0]
    Deltay = vec2[1]-vec1[1]
    r1 = m.sqrt(vec1[0]**2 + vec1[1]**2)
    r2 = m.sqrt(vec2[0]**2 + vec2[1]**2)
    dist = m.acosh(1 + 2*(Deltax**2 + Deltay**2)/((1 - r1**2)*(1-r2**2)))
    #print(vec1,vec2)
    #print(dist)
    return(m.acosh(1 + 2*(Deltax**2 + Deltay**2)/((1 - r1**2)*(1-r2**2))))

def draw_BTZ_geod(vec, rH): #takes two cartesian AdS points

    vec1 = vec[0]
    vec2 = vec[1]

    Nmax = 50000

    vec1p = CartesianToPolar(vec1)
    vec2p = CartesianToPolar(vec2)

    PoincareDiskPoints_C = (vec1, vec2)
    PoincareDiskPoints_P = (vec1p, vec2p)

    PoinGeodesicPoints = []

    if (vec2p[1] - vec1p[1]) % m.pi == 0: #geodesic connecting points is a line

        for n in range(0,Nmax):

            point = (vec1[0] + (n/Nmax)*(vec2[0] - vec1[0]), vec1[1] + (n/Nmax)*(vec2[1] - vec1[1]))

            PoinGeodesicPoints.append(point)      


    else: #geodesic connecting points is a circle

        PoincareDiskPoints = [] #Red and Blue points in polar coordinates

        angle_1 = min(PoincareDiskPoints_P[0][1], PoincareDiskPoints_P[1][1])
        angle_2 = max(PoincareDiskPoints_P[0][1], PoincareDiskPoints_P[1][1])
        angle = angle_2 - angle_1

        CircleData = PoincareGeodesicCircle(vec1, vec2)
                
        c_x = CircleData[0][0]
        c_y = CircleData[0][1]

        for n in range(0,Nmax): 

            Cxn = c_x + CircleData[1]*m.cos(2*m.pi*n/Nmax) 

            Cyn = c_y + CircleData[1]*m.sin(2*m.pi*n/Nmax) 

            if m.atan2(Cyn, Cxn) <= angle_2 and m.atan2(Cyn, Cxn) >= angle_1 and m.sqrt(Cxn**2 + Cyn**2) < 1.0: 

                PoinGeodesicPoints.append((Cxn,Cyn))

    BTZDiskGeodesicPoints=[]

    for vec in PoinGeodesicPoints: 

        BTZDiskGeodesicPoints.append(PolarToCartesian(PoincareDiskToBTZDisk(CartesianToPolar(vec),BTZDiskToBTZ(rH))))
                
    #returns the array of points defining the purple line (geodesic between vec1 and vec2 in BTZ.
    return(BTZDiskGeodesicPoints)

def draw_AdS_geod(vec):

    vec1 = vec[0]
    vec2 = vec[1]

    x_1 = vec1[0]
    x_2 = vec2[0]
    y_1 = vec1[1]
    y_2 = vec2[1]

    theta_ = m.atan2(y_2,x_2) - m.atan2(y_1,x_1)
	
    if theta_ > m.pi:

        X_1 = x_2
        Y_1 = y_2
        X_2 = x_1
        Y_2 = y_1
        #theta = 2*m.pi - theta_

    elif theta_ < -m.pi:

        X_1 = x_1
        Y_1 = y_1
        X_2 = x_2
        Y_2 = y_2
        #theta = 2*m.pi - abs(theta_)

    elif 0 < theta_ < m.pi:
	
        X_1 = x_1
        Y_1 = y_1
        X_2 = x_2
        Y_2 = y_2
        #theta = theta_

    else:

        X_1 = x_2
        Y_1 = y_2
        X_2 = x_1
        Y_2 = y_1
        #theta = abs(theta_)

    new_endpoints = ((X_1,Y_1),(X_2,Y_2))

    Circle_Data = PoincareGeodesicCircle(new_endpoints[0],new_endpoints[1])

    c_x = Circle_Data[0][0]
    c_y = Circle_Data[0][1]

    theta_2 = m.atan2(Y_2 - c_y, X_2 - c_x)*180/m.pi
    theta_1 = m.atan2(Y_1 - c_y, X_1 - c_x)*180/m.pi

    return(Circle_Data, theta_2, theta_1)

	
def CP(vec, rH): #takes and returns two cartestian AdS points

    vec1 = vec[0]
    vec2 = vec[1]

    vec1p = CartesianToPolar(vec1)
    vec2p = CartesianToPolar(vec2)

    PoincareDiskPoints_P = (vec1p, vec2p)

    BTZ_points = []

    for vec in PoincareDiskPoints_P:

        BTZ_points.append(PoincareDiskToBTZDisk(vec,BTZDiskToBTZ(rH)))

    if BTZ_points[0][1] < BTZ_points[1][1]:

        new_BTZ_points = ((BTZ_points[0][0], BTZ_points[0][1] + 2*m.pi), BTZ_points[1])

    else: #BTZ_points[0][1] >= BTZ_points[1][1]

        new_BTZ_points = (BTZ_points[0], (BTZ_points[1][0], BTZ_points[1][1] + 2*m.pi))

    new_AdS_points = []

    for vec in new_BTZ_points:

        new_AdS_points.append(PolarToCartesian(BTZDiskToPoincareDisk(vec, BTZDiskToBTZ(rH))))

    return(new_AdS_points)

cutoff = 0.0000001
max_r = 1 - cutoff

def EQ_phase(theta_A1, theta_A2, theta_B1, theta_B2, rH): #takes 4 polar angles in BTZ

    A_center = (theta_A1 + theta_A2)/2
    A_size = theta_A2 - theta_A1
    B_size = theta_B2 - theta_B1 + 2*m.pi

    vecA1 = (max_r, theta_A1)
    vecA2 = (max_r, theta_A2)
    vecB1 = (max_r, theta_B1)
    vecB2 = (max_r, theta_B2)

    BTZ_bdry_points_P = (vecA1, vecA2, vecB1, vecB2)	

    x_A1 = m.cos(theta_A1)*max_r
    y_A1 = m.sin(theta_A1)*max_r
    x_A2 = m.cos(theta_A2)*max_r
    y_A2 = m.sin(theta_A2)*max_r

    x_B1 = m.cos(theta_B1)*max_r
    y_B1 = m.sin(theta_B1)*max_r
    x_B2 = m.cos(theta_B2)*max_r
    y_B2 = m.sin(theta_B2)*max_r
		
    AdS_bdry_points_P = []

    for vec in BTZ_bdry_points_P:

        AdS_bdry_points_P.append(BTZDiskToPoincareDisk(vec, BTZDiskToBTZ(rH)))

    vecA1_AdS_P = AdS_bdry_points_P[0]
    vecA2_AdS_P = AdS_bdry_points_P[1]
    vecB1_AdS_P = AdS_bdry_points_P[2]
    vecB2_AdS_P = AdS_bdry_points_P[3]

    vecA1_AdS = PolarToCartesian(vecA1_AdS_P)
    vecA2_AdS = PolarToCartesian(vecA2_AdS_P)
    vecB1_AdS = PolarToCartesian(vecB1_AdS_P)
    vecB2_AdS = PolarToCartesian(vecB2_AdS_P)

    hor_max = m.sqrt((m.cosh(BTZDiskToBTZ(rH)*m.pi)-1)/(m.cosh(BTZDiskToBTZ(rH)*m.pi)+1))
    AdS_hpoint1 = (0,-hor_max)
    AdS_hpoint2 = (0,hor_max)

    horizon_area = AdS_distance((AdS_hpoint1, AdS_hpoint2))

    c_1 = AdS_distance((vecA2_AdS, vecB1_AdS)) + AdS_distance((vecB2_AdS, vecA1_AdS)) + horizon_area
    c_2 = AdS_distance((vecA1_AdS, vecA2_AdS)) + AdS_distance(CP((vecB1_AdS, vecB2_AdS), rH))
    c_3 = AdS_distance((vecA2_AdS, vecB1_AdS)) + AdS_distance(CP((vecB2_AdS, vecA1_AdS), rH))
    c_4 = AdS_distance((vecB2_AdS, vecA1_AdS)) + AdS_distance(CP((vecB1_AdS, vecA2_AdS), rH))

    #print(c_1,c_2,c_3,c_4)

    EW_type = np.argmin((c_1, c_2, c_3, c_4))

    if EW_type == 2 or EW_type == 3:

        return 'BH outside EW'

    elif EW_type == 1:

        return 'EW disconn'

    #Define objective function to be minimized
	
    def f_CM(opt_points): #two angles in radians, and two y-coordinates for horizon points

        geod_point1 = (geod1_data[0][0][0] + geod1_data[0][1]*m.cos(opt_points[0]), geod1_data[0][0][1] + geod1_data[0][1]*m.sin(opt_points[0]))
        geod_point2 = (geod2_data[0][0][0] + geod2_data[0][1]*m.cos(opt_points[1]), geod2_data[0][0][1] + geod2_data[0][1]*m.sin(opt_points[1]))

        horizon_point1 = (0, opt_points[2])
        horizon_point2 = (0, opt_points[3])

        #S_A--------------------------------------------------------------

        S_A = min(AdS_distance((vecA1_AdS, vecA2_AdS)), AdS_distance(CP((vecA1_AdS, vecA2_AdS), rH)) + horizon_area)

        #S_B--------------------------------------------------------------

        S_B = min(AdS_distance(CP((vecB1_AdS, vecB2_AdS), rH)), AdS_distance((vecB1_AdS, vecB2_AdS)) + horizon_area)

		#S_Aa-------------------------------------------------------------

        if horizon_point1 == horizon_point2:

            if horizon_point1[1] < 0: #horizon belongs to a

                S_Aa = min(AdS_distance((geod_point1, geod_point2)) + horizon_area, AdS_distance(CP((geod_point1, geod_point2), rH)))

            else: # horizon_point1[1] >= 0: #horizon belongs to b

                S_Aa = min(AdS_distance((geod_point1, geod_point2)), AdS_distance(CP((geod_point1, geod_point2), rH)) + horizon_area)

        else:

            if horizon_point1[1] < horizon_point2[1]: #horizon part of a DOES NOT cross critical point

                S_Aa = min(AdS_distance((geod_point1, horizon_point1)) + AdS_distance((geod_point2, horizon_point2)), AdS_distance((horizon_point1, horizon_point2)) + AdS_distance((geod_point1, geod_point2)), AdS_distance(CP((horizon_point1, horizon_point2), rH)) + AdS_distance(CP((geod_point1, geod_point2), rH)))

            else: # horizon_point1[1] > horizon_point2[1]: #horizon part of a crosses critical point

                S_Aa = min(AdS_distance((horizon_point1, horizon_point2)) + AdS_distance(CP((geod_point1, geod_point2), rH)), AdS_distance(CP((horizon_point1, horizon_point2), rH)) + AdS_distance((geod_point1, geod_point2)))

		#S_Ba-------------------------------------------------------------

        if horizon_point1 == horizon_point2:

            if horizon_point1[1] < 0: #horizon belongs to a

                S_Ba = min(AdS_distance((vecA1_AdS, vecA2_AdS)) + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS)), AdS_distance(CP((vecA1_AdS, vecA2_AdS), rH)) + horizon_area + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS)), AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance(CP((vecB1_AdS, vecB2_AdS), rH)) + horizon_area, AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance((vecB1_AdS, vecB2_AdS)))

            else: # horizon_point1[1] >= 0: #horizon belongs to b

                S_Ba = min(AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance(CP((vecB1_AdS, vecB2_AdS), rH)), AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance((vecB1_AdS, vecB2_AdS)) + horizon_area, AdS_distance((vecA1_AdS, vecA2_AdS)) + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS)) + horizon_area, AdS_distance(CP((vecA1_AdS, vecA2_AdS), rH)) + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS)))

        else:

            if horizon_point2[1] > horizon_point1[1]: #horizon part of a DOES NOT cross critical point

                S_Ba = min(AdS_distance((vecA1_AdS, vecA2_AdS)) + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS)) + AdS_distance(CP((horizon_point1, horizon_point2), rH)), AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance(CP((vecB1_AdS, vecB2_AdS), rH)) + AdS_distance((horizon_point1, horizon_point2)), AdS_distance(CP((vecA1_AdS, vecA2_AdS), rH)) + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS)) + AdS_distance((horizon_point1, horizon_point2)), AdS_distance((vecB1_AdS, vecB2_AdS)) + AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance(CP((horizon_point1, horizon_point2), rH)))

            else: # horizon_point2[1] < horizon_point1[1]: #horizon part of a crosses critical point

                S_Ba = min(AdS_distance((vecA1_AdS, vecA2_AdS)) + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS)) + AdS_distance((horizon_point1, horizon_point2)), AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance(CP((vecB1_AdS, vecB2_AdS), rH)) + AdS_distance(CP((horizon_point1, horizon_point2), rH)), AdS_distance(CP((vecA1_AdS, vecA2_AdS), rH)) + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS)) + AdS_distance(CP((horizon_point1, horizon_point2), rH)), AdS_distance((vecB1_AdS, vecB2_AdS)) + AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance((horizon_point1, horizon_point2)))
	          
        f_CM = S_A + S_B + S_Aa - S_Ba

        return f_CM

    #minimize the objective

    geod1_data = draw_AdS_geod((vecA1_AdS, vecB2_AdS))
    geod2_data = draw_AdS_geod((vecA2_AdS, vecB1_AdS))

    lb_1 = (geod1_data[1]*m.pi/180)%(2*m.pi)
    ub_1 = (geod1_data[2]*m.pi/180)%(2*m.pi)

    lb_2 = (geod2_data[1]*m.pi/180)%(2*m.pi)
    ub_2 = (geod2_data[2]*m.pi/180)%(2*m.pi)
    
    lb_3 = -hor_max
    ub_3 = hor_max

    bnds = opt.Bounds([lb_1, lb_2, lb_3, lb_3], [ub_1, ub_2, ub_3 - 0.0000001, ub_3 - 0.0000001], keep_feasible = True)
    constraint_1 = opt.LinearConstraint([[0,0,1,-1],[0,0,1,0]], [0,lb_3],[0,0], keep_feasible = True)
    constraint_2 = opt.LinearConstraint([[0,0,1,-1],[0,0,1,0]], [0,0],[0,ub_3], keep_feasible = True)
    constraint_3 = opt.LinearConstraint([[0,0,-1,1]], [0.000000001],[2*ub_3], keep_feasible = True)
    constraint_4 = opt.LinearConstraint([[0,0,1,-1]], [0.000000001],[2*ub_3], keep_feasible = True)

    E_CM = None

    for i in range(20):

        A = [rand.randint(1, 999), rand.randint(1, 999), rand.randint(1, 999), rand.randint(1, 999)]

        theta_1 = lb_1 + (A[0]/1000)*(ub_1 - lb_1)
        theta_2 = lb_2 + (A[1]/1000)*(ub_2 - lb_2)

###################Constraint 1#############################     
        if A[2] > 500:
            P = 1000- A[2]
        else:
            P = A[2]
            
        h_1 = lb_3 + (P/1000)*(ub_3 - lb_3)
        h_2 = lb_3 + (P/1000)*(ub_3 - lb_3)
        
        f__CM = opt.minimize(f_CM, [theta_1, theta_2, h_1, h_2], bounds = bnds, constraints = constraint_1)
    
        if i == 0:

            E_CM = f__CM
            domain = 1

        else:
    
            if f__CM.fun < E_CM.fun:
                    
                E_CM = f__CM
                domain = 1

###################Constraint 2#############################                      

        if A[2] < 500:
            
            Q = 1000 - A[2]
            
        else:
            
            Q = A[2]
            
        h_1 = lb_3 + (Q/1000)*(ub_3 - lb_3)
        h_2 = lb_3 + (Q/1000)*(ub_3 - lb_3)

        f__CM = opt.minimize(f_CM, [theta_1, theta_2, h_1, h_2], bounds = bnds, constraints = constraint_2)
    
        if f__CM.fun < E_CM.fun:
                
            E_CM = f__CM
            domain = 2

###################Constraint 3#############################     


        if A[2] == A[3]:

            A[3] = (A[3]- rand.randint(1, 999))%1000

        if A[3] < A[2]:

            R = A[2]
            S = A[3]
            
        else:

            R = A[3]
            S = A[2]
            
            
        h_1 = lb_3 + (S/1000)*(ub_3 - lb_3)
        h_2 = lb_3 + (R/1000)*(ub_3 - lb_3)

        f__CM = opt.minimize(f_CM, [theta_1, theta_2, h_1, h_2], bounds = bnds, constraints = constraint_3)

        if f__CM.fun < E_CM.fun:
                
            E_CM = f__CM
            domain = 3

###################Constraint 4#############################                 

        if A[3] < A[2]:

            T = A[2]
            U = A[3]
            
        else:

            T = A[3]
            U = A[2]
            
            
        h_1 = lb_3 + (T/1000)*(ub_3 - lb_3)
        h_2 = lb_3 + (U/1000)*(ub_3 - lb_3)

        f__CM = opt.minimize(f_CM, [theta_1, theta_2, h_1, h_2], bounds = bnds, constraints = constraint_4)

        if f__CM.fun < E_CM.fun:
                
            E_CM = f__CM
            domain = 4


    

    #print(E_CM.x[2],E_CM.x[3])
    #print(E_CM.fun)
    
    A = f_CM((E_CM.x[0], E_CM.x[1], E_CM.x[2], E_CM.x[2]))
    B = f_CM((E_CM.x[0], E_CM.x[1], -E_CM.x[2], -E_CM.x[2]))

    #print(A)
    #print(B)

    if A <= E_CM.fun:

        E_CM.x[3] = E_CM.x[2]
        E_CM.fun = A

    if B <= E_CM.fun:

        E_CM.x[2] = -E_CM.x[2]
        E_CM.x[3] = E_CM.x[2]
        E_CM.fun = B

    #print(E_CM.x[2],E_CM.x[3])

    geod_point1 =(geod1_data[0][0][0] + geod1_data[0][1]*m.cos(E_CM.x[0]), geod1_data[0][0][1] + geod1_data[0][1]*m.sin(E_CM.x[0]))
    geod_point2 =(geod2_data[0][0][0] + geod2_data[0][1]*m.cos(E_CM.x[1]), geod2_data[0][0][1] + geod2_data[0][1]*m.sin(E_CM.x[1]))
    horizon_point1 = (0, E_CM.x[2])
    horizon_point2 = (0, E_CM.x[3])

    min_points_C = [geod_point1, geod_point2, horizon_point1, horizon_point2]

    if E_CM.x[2] == E_CM.x[3]:
        

        if E_CM.x[2] < 0:
            #print(1)
            cand = np.array([AdS_distance((vecA1_AdS, vecA2_AdS)) + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS)), AdS_distance(CP((vecA1_AdS, vecA2_AdS), rH)) + horizon_area + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS)), AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance(CP((vecB1_AdS, vecB2_AdS), rH)) + horizon_area, AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance((vecB1_AdS, vecB2_AdS))])

            #print(cand)

            sorted_cand = np.argsort(cand)

            minima = (sorted_cand[0],sorted_cand[1])

        else:
            #print(2)
            cand = np.array([AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance(CP((vecB1_AdS, vecB2_AdS), rH)), AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance((vecB1_AdS, vecB2_AdS)) + horizon_area, AdS_distance((vecA1_AdS, vecA2_AdS)) + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS)) + horizon_area, AdS_distance(CP((vecA1_AdS, vecA2_AdS), rH)) + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS))])

            #print(cand)

            sorted_cand = np.argsort(cand)

            minima = (sorted_cand[0],sorted_cand[1])

        if abs(cand[minima[0]] - cand[minima[1]]) < 0.00001:

            if minima == (1,2) or minima == (2,1):

                phase = 1

            if minima == (0,2) or minima == (2,0):

                phase = 2

            if minima == (0,3) or minima == (3,0):

                phase = 3

    elif E_CM.x[2] < E_CM.x[3]:

        #print(3)

        Aa_cand = np.argmin((AdS_distance((geod_point1, horizon_point1)) + AdS_distance((geod_point2, horizon_point2)), AdS_distance((horizon_point1, horizon_point2)) + AdS_distance((geod_point1, geod_point2)), AdS_distance(CP((horizon_point1, horizon_point2), rH)) + AdS_distance(CP((geod_point1, geod_point2), rH))))

        if Aa_cand != 0:

            if abs(E_CM.x[2] + E_CM.x[3]) > 0.0000001:

                phase = 3

            else:

                phase = 2

        else:

            phase = 4

    elif E_CM.x[2] > E_CM.x[3]:
        #print(4)

        if abs(E_CM.x[2] + E_CM.x[3]) > 0.0000001:

            phase = 3

        else:

            phase = 2


    #else:

        #phase = 4

    #cand1 = np.array([AdS_distance((vecA1_AdS, vecA2_AdS)) + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS)), AdS_distance(CP((vecA1_AdS, vecA2_AdS), rH)) + horizon_area + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS)), AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance(CP((vecB1_AdS, vecB2_AdS), rH)) + horizon_area, AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance((vecB1_AdS, vecB2_AdS))])

    #print(cand1)

    #cand2 = np.array([AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance(CP((vecB1_AdS, vecB2_AdS), rH)), AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance((vecB1_AdS, vecB2_AdS)) + horizon_area, AdS_distance((vecA1_AdS, vecA2_AdS)) + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS)) + horizon_area, AdS_distance(CP((vecA1_AdS, vecA2_AdS), rH)) + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS))])

    #print(cand2)

    #sorted_cand1 = np.argsort(cand1)

    #sorted_cand2 = np.argsort(cand2)

    #minima1 = (sorted_cand1[0],sorted_cand1[1])

    #minima2 = (sorted_cand2[0],sorted_cand2[1])

    #if abs(cand1[minima1[0]] - cand1[minima1[1]]) < 0.0001 and abs(cand2[minima2[0]] - cand2[minima2[1]]) > 0.0001:

        #if minima1 == (1,2) or minima1 == (2,1):

            #phase = 1

        #if minima1 == (0,2) or minima1 == (2,0):

            #phase = 2

        #if minima1 == (0,3) or minima1 == (3,0):

            #phase = 3

    #elif abs(cand2[minima2[0]] - cand2[minima2[1]]) < 0.0001 and abs(cand1[minima1[0]] - cand1[minima1[1]]) > 0.0001:

        #if minima2 == (1,2) or minima2 == (2,1):

            #phase = 1

        #if minima2 == (0,2) or minima2 == (2,0):

            #phase = 2

        #if minima2 == (0,3) or minima2 == (3,0):

            #phase = 3

    #if abs(cand1[minima1[0]] - cand1[minima1[1]]) < 0.0001:

        #if minima1 == (1,2) or minima1 == (2,1):

            #phase = 1

        #if minima1 == (0,2) or minima1 == (2,0):

            #phase = 2

        #if minima1 == (0,3) or minima1 == (3,0):

            #phase = 3

    #elif abs(cand2[minima2[0]] - cand2[minima2[1]]) < 0.0001:

        #if minima2 == (1,2) or minima2 == (2,1):

            #phase = 1

        #if minima2 == (0,2) or minima2 == (2,0):

            #phase = 2

        #if minima2 == (0,3) or minima2 == (3,0):

            #phase = 3

    #else:

        #phase = 4

    

    return phase

#data = np.array([])
#i = 0

#for j in range(100):

    #A = EQ_phase(-12*m.pi/40,12*m.pi/40,14*m.pi/40,-14*m.pi/40,0.28)

    #np.append(data, int(A[0]))

    #if A[0] == 3:

        #i = i+1

    #if A[0] == 3:
        
        #break

#print(A[1])

#print(i)
        





#print(EQ_phase(-1*m.pi/12, 3*m.pi/12, 4*m.pi/12, -4*m.pi/12, 0.3))


file = open('dict_EQ','rb')
Dict = pickle.load(file)
file.close()

new_dict_entries = {}

rH = 0.19

A_center = 0

for A_size in range(1,40):

    for B_size in range(1,40):

        theta_A1 = (A_center - A_size)*m.pi/40
        theta_A2 = (A_center + A_size)*m.pi/40
        theta_B1 = (40 - B_size)*m.pi/40
        theta_B2 = (-40 + B_size)*m.pi/40

        if theta_A2 < theta_B1:
            
            A = EQ_phase(theta_A1, theta_A2, theta_B1, theta_B2, rH)

            if Dict[(A_center, A_size, B_size, rH)] != A:

                votes = []

                for i in range(15):

                    vote = EQ_phase(theta_A1, theta_A2, theta_B1, theta_B2, rH)
                    print(vote)
                    votes.append(vote)

                counts = np.bincount(votes)
                new_dict_entries[(A_center, A_size, B_size, rH)] = np.argmax(counts)

for x in new_dict_entries:

    Dict[x] = new_dict_entries[x]

file = open('dict_EQ','wb')
pickle.dump(Dict, file)
file.close()
                
                    

#new_dict_entries = {}

#for i in range(1,40):

#    A_center = i

#    A_size = A_center + 0.0001

#    for B_size in range(1,40):

#        theta_A1 = (A_center - A_size)*m.pi/40
#        theta_A2 = (A_center + A_size)*m.pi/40
#        theta_B1 = (40 - B_size)*m.pi/40
#        theta_B2 = (-40 + B_size)*m.pi/40

#        if theta_A2 < theta_B1:

#            new_dict_entries[(A_center, int(A_size), B_size, rH)] = EQ_phase(theta_A1, theta_A2, theta_B1, theta_B2, rH)


#file = open('dict_EQ','rb')
#Dict = pickle.load(file)
#file.close()

#for x in new_dict_entries:

#    Dict[x] = new_dict_entries[x]

#file = open('dict_EQ','wb')
#pickle.dump(Dict, file)
#file.close()



            





