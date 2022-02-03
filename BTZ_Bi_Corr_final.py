#!/usr/bin/env python3
import numpy as np
import math as m
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as pat
import scipy.optimize as opt
import random as rand
from scipy.interpolate import interp1d
from scipy import interpolate
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

    BTZDiskGeodesicPoints_x=[]
    BTZDiskGeodesicPoints_y=[]

    for vec in PoinGeodesicPoints: 

        BTZDiskGeodesicPoints_x.append(PolarToCartesian(PoincareDiskToBTZDisk(CartesianToPolar(vec),BTZDiskToBTZ(rH)))[0])
        BTZDiskGeodesicPoints_y.append(PolarToCartesian(PoincareDiskToBTZDisk(CartesianToPolar(vec),BTZDiskToBTZ(rH)))[1])
                
    #returns the array of points defining the purple line (geodesic between vec1 and vec2 in BTZ.

    tck, u = interpolate.splprep([BTZDiskGeodesicPoints_x,BTZDiskGeodesicPoints_y],s=0)
    unew = np.arange(0, 1, 1/50000)
    out = interpolate.splev(unew, tck)
   
    #f = interp1d(BTZDiskGeodesicPoints_x, BTZDiskGeodesicPoints_y,kind = 'cubic')
      
    return(out)

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

def Bi_Corr_Meas(theta_A1, theta_A2, theta_B1, theta_B2, rH): #takes 4 polar angles in BTZ

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

    fig, ax = plt.subplots()

    ax.set_xlim((-1.1,1.1))
    ax.set_ylim((-1.1,1.1))

    ax.set_aspect('equal')
	
    ax.tick_params(top = 'off', bottom = 'off', left = 'off', right = 'off', labeltop = 'off', labelbottom = 'off', labelleft = 'off', labelright = 'off')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    #draw boundary	

    boundary = plt.Circle((0,0), 1, color = 'black', fill = False, clip_on = False)
    ax.add_artist(boundary)

    #draw horizon

    hor_max = m.sqrt((m.cosh(BTZDiskToBTZ(rH)*m.pi)-1)/(m.cosh(BTZDiskToBTZ(rH)*m.pi)+1))
    AdS_hpoint1 = (0,-hor_max)
    AdS_hpoint2 = (0,hor_max)
    
    horizon = plt.Circle((0,0), rH, color = 'black', fill = False, clip_on = False)
    ax.add_artist(horizon)

    #plot boundary regions in BTZ		

    arc1 = pat.Arc((0,0), 2, 2, angle = 0, theta1 = 180*theta_A1/m.pi, theta2 = 180*theta_A2/m.pi, linewidth = 2, color = 'blue')
    arc2 = pat.Arc((0,0), 2, 2, angle = 0, theta1 = 180*theta_B1/m.pi, theta2 = 180*theta_B2/m.pi, linewidth = 2, color = 'blue')
    ax.add_artist(arc1)
    #ax.add_artist(arc2)

    #plot boundary regions in AdS
	
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

    #draw entanglement wedge and define I_AB

    horizon_area = AdS_distance((AdS_hpoint1, AdS_hpoint2))

    c_1 = AdS_distance((vecA2_AdS, vecB1_AdS)) + AdS_distance((vecB2_AdS, vecA1_AdS)) + horizon_area
    c_2 = AdS_distance((vecA1_AdS, vecA2_AdS)) + AdS_distance(CP((vecB1_AdS, vecB2_AdS), rH))
    c_3 = AdS_distance((vecA2_AdS, vecB1_AdS)) + AdS_distance(CP((vecB2_AdS, vecA1_AdS), rH))
    c_4 = AdS_distance((vecB2_AdS, vecA1_AdS)) + AdS_distance(CP((vecB1_AdS, vecA2_AdS), rH))

    EW_type = np.argmin((c_1, c_2, c_3, c_4))

    geod1_data = draw_AdS_geod((vecA1_AdS, vecB2_AdS))
    geod2_data = draw_AdS_geod((vecA2_AdS, vecB1_AdS))

    if EW_type == 0:
	
        #axs[SP_BTZ[0],SP_BTZ[1]].scatter(*zip(*draw_BTZ_geod((vecA1_AdS, vecB2_AdS), rH)), s=1, color = 'purple') 
        #axs[SP_BTZ[0],SP_BTZ[1]].scatter(*zip(*draw_BTZ_geod((vecA2_AdS, vecB1_AdS), rH)), s=1, color = 'purple')
        #axs[SP_BTZ[0],SP_BTZ[1]].scatter(*zip(*draw_BTZ_geod((AdS_hpoint1, AdS_hpoint2), rH)), s=1, color = 'purple')

        #axs[SP_BTZ[0],SP_BTZ[1]].plot(draw_BTZ_geod((vecA1_AdS, vecB2_AdS), rH)[0],draw_BTZ_geod((vecA1_AdS, vecB2_AdS), rH)[1], color = 'purple', linewidth = 2)
        #axs[SP_BTZ[0],SP_BTZ[1]].plot(draw_BTZ_geod((vecA2_AdS, vecB1_AdS), rH)[0],draw_BTZ_geod((vecA2_AdS, vecB1_AdS), rH)[1], color = 'purple', linewidth = 2)
        #axs[SP_BTZ[0],SP_BTZ[1]].plot(draw_BTZ_geod((AdS_hpoint1, AdS_hpoint2), rH)[0],draw_BTZ_geod((AdS_hpoint1, AdS_hpoint2), rH)[1], color = 'black', linewidth = 0.8)

        AdS_geod1 = pat.Arc((geod1_data[0][0][0],geod1_data[0][0][1]), 2*geod1_data[0][1], 2*geod1_data[0][1], angle = 0, theta1 = geod1_data[1], theta2 = geod1_data[2], linewidth = 3, color = 'purple')
        AdS_geod2 = pat.Arc((geod2_data[0][0][0],geod2_data[0][0][1]), 2*geod2_data[0][1], 2*geod2_data[0][1], angle = 0, theta1 = geod2_data[1], theta2 = geod2_data[2], linewidth = 3, color = 'purple')
        
    elif EW_type == 1:

        print('Entanglement wedge disconnected')
        #return()

    elif EW_type == 2 or EW_type == 3:

        print('Entanglement wedge connected but does not contain black hole')
        #return()

    S_A = min(AdS_distance((vecA1_AdS, vecA2_AdS)), AdS_distance(CP((vecA1_AdS, vecA2_AdS), rH)) + horizon_area)
    S_B = min(AdS_distance(CP((vecB1_AdS, vecB2_AdS), rH)), AdS_distance((vecB1_AdS, vecB2_AdS)) + horizon_area)
    S_AB = min(c_1,c_2,c_3,c_4)

    I_AB = S_A + S_B - S_AB

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

        #S_a--------------------------------------------------------------

        if horizon_point1 == horizon_point2:

            if horizon_point1[1] < 0: #horizon belongs to a

                S_a = AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + horizon_area

            else: #horizon_point1[1] >= 0: #horizon belongs to b

                S_a = AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2))

        else:

            if horizon_point2[1] > horizon_point1[1]: #horizon part of a DOES NOT cross critical point

                S_a = AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance((horizon_point1, horizon_point2))

            else: # horizon_point2[1] < horizon_point1[1]: #horizon part of a crosses critical point

                S_a = AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance(CP((horizon_point2, horizon_point1), rH))

		#S_AB-------------------------------------------------------------

        S_AB = min(c_1,c_2,c_3,c_4)

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

                S_Ba = min(AdS_distance((vecA1_AdS, vecA2_AdS)) + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS)) + horizon_area, AdS_distance(CP((vecA1_AdS, vecA2_AdS), rH)) + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS)), AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance(CP((vecB1_AdS, vecB2_AdS), rH)), AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance((vecB1_AdS, vecB2_AdS)) + horizon_area)

        else:

            if horizon_point2[1] > horizon_point1[1]: #horizon part of a DOES NOT cross critical point

                S_Ba = min(AdS_distance((vecA1_AdS, vecA2_AdS)) + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS)) + AdS_distance(CP((horizon_point1, horizon_point2), rH)), AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance(CP((vecB1_AdS, vecB2_AdS), rH)) + AdS_distance((horizon_point1, horizon_point2)), AdS_distance(CP((vecA1_AdS, vecA2_AdS), rH)) + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS)) + AdS_distance((horizon_point1, horizon_point2)), AdS_distance((vecB1_AdS, vecB2_AdS)) + AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance(CP((horizon_point1, horizon_point2), rH)))

            else: # horizon_point2[1] < horizon_point1[1]: #horizon part of a crosses critical point

                S_Ba = min(AdS_distance((vecA1_AdS, vecA2_AdS)) + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS)) + AdS_distance((horizon_point1, horizon_point2)), AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance(CP((vecB1_AdS, vecB2_AdS), rH)) + AdS_distance(CP((horizon_point1, horizon_point2), rH)), AdS_distance(CP((vecA1_AdS, vecA2_AdS), rH)) + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS)) + AdS_distance(CP((horizon_point1, horizon_point2), rH)), AdS_distance((vecB1_AdS, vecB2_AdS)) + AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance((horizon_point1, horizon_point2)))

		#S_ABa------------------------------------------------------------

        if horizon_point1 == horizon_point2:

            if horizon_point1[1] < 0: #horizon belongs to a

                S_ABa = AdS_distance((vecB1_AdS, geod_point2)) + AdS_distance((vecB2_AdS, geod_point1))

            else: # horizon_point1[1] >= 0: #horizon belongs to b

                S_ABa = AdS_distance((vecB1_AdS, geod_point2)) + AdS_distance((vecB2_AdS, geod_point1)) + horizon_area

        else:

            if horizon_point2[1] > horizon_point1[1]: #horizon part of b crosses critical point

                S_ABa = AdS_distance((vecB1_AdS, geod_point2)) + AdS_distance((vecB2_AdS, geod_point1)) + AdS_distance(CP((horizon_point1, horizon_point2), rH))

            else: # horizon_point2[1] < horizon_point1[1]: #horizon part of b DOES NOT cross critical point

                S_ABa = AdS_distance((vecB1_AdS, geod_point2)) + AdS_distance((vecB2_AdS, geod_point1)) + AdS_distance((horizon_point1, horizon_point2))

	
        S_Bab = S_A

        S_Aab = S_B

        S_ABb = S_a

        S_ab = S_AB

        S_Bb = S_Aa

        S_Ab = S_Ba

        S_b = S_ABa

	

        #I_AB-------------------------------------------

        #f_CM = S_A + S_B - S_AB

        #E_sq-------------------------------------------

        #f_CM = S_Aa + S_Ba - S_ABa - S_a

        #E_P--------------------------------------------

        #f_CM = S_Aa

        #E_R--------------------------------------------

        #f_CM = S_AB + 2*S_Aa - S_a - S_ABa
		
        #E_Q--------------------------------------------

        f_CM = S_A + S_B + S_Aa - S_Ba

        #f_CM = S_Aa

        return f_CM

    #minimize the objective

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
    
    #-----------------------------------
    #------THIS PART FOR FIGURES--------
    #-----------------------------------   
    #E_CM.x[0] = -2.3*m.pi/2
    #E_CM.x[1] = 2.3*m.pi/2
    #E_CM.x[2] = -hor_max/2 
    #E_CM.x[3] = hor_max/2
    #-----------------------------------
    #-----------------------------------
    #-----------------------------------

    geod_point1 =(geod1_data[0][0][0] + geod1_data[0][1]*m.cos(E_CM.x[0]), geod1_data[0][0][1] + geod1_data[0][1]*m.sin(E_CM.x[0]))
    geod_point2 =(geod2_data[0][0][0] + geod2_data[0][1]*m.cos(E_CM.x[1]), geod2_data[0][0][1] + geod2_data[0][1]*m.sin(E_CM.x[1]))
    horizon_point1 = (0, E_CM.x[2])
    horizon_point2 = (0, E_CM.x[3])

    gp1_AdS_P = CartesianToPolar(geod_point1)
    gp2_AdS_P = CartesianToPolar(geod_point2)

    gp1_BTZ_P = PoincareDiskToBTZDisk(gp1_AdS_P,BTZDiskToBTZ(rH))
    gp2_BTZ_P = PoincareDiskToBTZDisk(gp2_AdS_P,BTZDiskToBTZ(rH))

    gp1_BTZ = PolarToCartesian(gp1_BTZ_P)
    gp2_BTZ = PolarToCartesian(gp2_BTZ_P)

    hp1_AdS_P = CartesianToPolar(horizon_point1)
    hp2_AdS_P = CartesianToPolar(horizon_point2)

    hp1_BTZ_P = PoincareDiskToBTZDisk(hp1_AdS_P,BTZDiskToBTZ(rH))
    hp2_BTZ_P = PoincareDiskToBTZDisk(hp2_AdS_P,BTZDiskToBTZ(rH))

    hp1_BTZ = PolarToCartesian(hp1_BTZ_P)
    hp2_BTZ = PolarToCartesian(hp2_BTZ_P)

    min_points_C = [geod_point1, geod_point2, horizon_point1, horizon_point2]
    

    def surf_CM(x):

        geod_point1 = (geod1_data[0][0][0] + geod1_data[0][1]*m.cos(x[0]), geod1_data[0][0][1] + geod1_data[0][1]*m.sin(x[0]))
        geod_point2 = (geod2_data[0][0][0] + geod2_data[0][1]*m.cos(x[1]), geod2_data[0][0][1] + geod2_data[0][1]*m.sin(x[1]))
        horizon_point1 = (0, x[2])
        horizon_point2 = (0, x[3])        

        v0 = np.array([1,0,0,0,0,0,0,0,0,0,0,0,0,0])#A1-A2
        v1 = np.array([0,1,0,0,0,0,0,0,0,0,0,0,0,0])#A1-A2 CP
        v2 = np.array([0,0,1,0,0,0,0,0,0,0,0,0,0,0])#B1-B2
        v3 = np.array([0,0,0,1,0,0,0,0,0,0,0,0,0,0])#B1-B2 CP
        v4 = np.array([0,0,0,0,1,0,0,0,0,0,0,0,0,0])#A2-G2
        v5 = np.array([0,0,0,0,0,1,0,0,0,0,0,0,0,0])#B1-G2
        v6 = np.array([0,0,0,0,0,0,1,0,0,0,0,0,0,0])#A1-G1
        v7 = np.array([0,0,0,0,0,0,0,1,0,0,0,0,0,0])#B2-G1
        v8 = np.array([0,0,0,0,0,0,0,0,1,0,0,0,0,0])#G1-G2
        v9 = np.array([0,0,0,0,0,0,0,0,0,1,0,0,0,0])#G1-G2 CP
        v10 = np.array([0,0,0,0,0,0,0,0,0,0,1,0,0,0])#G1-H1
        v11 = np.array([0,0,0,0,0,0,0,0,0,0,0,1,0,0])#G2-H2
        v12 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,1,0])#H1-H2
        v13 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,1])#H1-H2 CP

        #A--------------------------------------------------------------

        surf_A = np.argmin((AdS_distance((vecA1_AdS, vecA2_AdS)), AdS_distance(CP((vecA1_AdS, vecA2_AdS), rH)) + horizon_area))

        if surf_A == 0:

            A = v0

        else:

            A = v1 + v12 + v13

        #B--------------------------------------------------------------

        surf_B = np.argmin((AdS_distance(CP((vecB1_AdS, vecB2_AdS), rH)), AdS_distance((vecB1_AdS, vecB2_AdS)) + horizon_area))

        if surf_B == 0:

            B = v3

        else:

            B = v2 + v12 + v13

        #a--------------------------------------------------------------

        if horizon_point1 == horizon_point2:

            if horizon_point1[1] < 0: #horizon belongs to a

                a = v4 + v6 + v12 + v13

            else: # horizon_point1[1] >= 0: #horizon belongs to b

                a = v4 + v6

        else:

            if horizon_point2[1] > horizon_point1[1]: #horizon part of a DOES NOT cross critical point

                a = v4 + v6 + v12

            else: # horizon_point2[1] < horizon_point1[1]: #horizon part of a crosses critical point

                a = v4 + v6 + v13

        #AB-------------------------------------------------------------

        AB = v6 + v7 + v4 + v5 + v12 + v13

        #Aa-------------------------------------------------------------

        if horizon_point1 == horizon_point2:

            if horizon_point1[1] < 0: #horizon belongs to a

                surf_Aa = np.argmin((AdS_distance((geod_point1, geod_point2)) + horizon_area, AdS_distance(CP((geod_point1, geod_point2), rH))))

                if surf_Aa == 0:

                    Aa = v8 + v12 + v13

                else:

                    Aa = v9

            else: # horizon_point1[1] >= 0: #horizon belongs to b

                    
                surf_Aa = np.argmin((AdS_distance((geod_point1, geod_point2)), AdS_distance(CP((geod_point1, geod_point2), rH)) + horizon_area))

                if surf_Aa == 0:

                    Aa = v8

                else:

                    Aa = v9 + v12 + v13

        else:

            if horizon_point1[1] < horizon_point2[1]: #horizon part of a DOES NOT cross critical point

                surf_Aa = np.argmin((AdS_distance((geod_point1, horizon_point1)) + AdS_distance((geod_point2, horizon_point2)), AdS_distance((horizon_point1, horizon_point2)) + AdS_distance((geod_point1, geod_point2)), AdS_distance(CP((horizon_point1, horizon_point2), rH)) + AdS_distance(CP((geod_point1, geod_point2), rH))))

                if surf_Aa == 0:

                    Aa = v10 + v11

                elif surf_Aa == 1:

                    Aa = v8 + v12

                else: # surf_Aa == 2:

                    Aa = v9 + v13

            else: # horizon_point1[1] > horizon_point2[1]: #horizon part of a crosses critical point

                surf_Aa = np.argmin((AdS_distance((horizon_point1, horizon_point2)) + AdS_distance(CP((geod_point1, geod_point2), rH)), AdS_distance(CP((horizon_point1, horizon_point2), rH)) + AdS_distance((geod_point1, geod_point2))))

                if surf_Aa == 0:

                    Aa = v9 + v12

                else:

                    Aa = v8 + v13

            #Ba-------------------------------------------------------------



        if horizon_point1 == horizon_point2:

            if horizon_point1[1] < 0: #horizon belongs to a

                surf_Ba = np.argmin((AdS_distance((vecA1_AdS, vecA2_AdS)) + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS)), AdS_distance(CP((vecA1_AdS, vecA2_AdS), rH)) + horizon_area + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS)), AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance(CP((vecB1_AdS, vecB2_AdS), rH)) + horizon_area, AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance((vecB1_AdS, vecB2_AdS))))
                if surf_Ba == 0:

                    Ba = v0 + v5 + v7

                elif surf_Ba == 1:

                    Ba = v1 + v12 + v13 + v5 + v7

                elif surf_Ba == 2:

                    Ba = v4 + v6 + v3 + v12 + v13

                else:

                    Ba = v4 + v6 + v2

            else: # horizon_point1[1] >= 0: #horizon belongs to b

                surf_Ba = np.argmin((AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance(CP((vecB1_AdS, vecB2_AdS), rH)), AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance((vecB1_AdS, vecB2_AdS)) + horizon_area, AdS_distance((vecA1_AdS, vecA2_AdS)) + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS)) + horizon_area, AdS_distance(CP((vecA1_AdS, vecA2_AdS), rH)) + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS))))

                if surf_Ba == 0:

                    Ba = v4 + v6 + v3

                elif surf_Ba == 1:

                    Ba = v4 + v6 + v2 + v12 + v13

                elif surf_Ba == 2:

                    Ba = v0 + v5 + v7 + v12 + v13

                else:

                    Ba = v1 + v5 + v7

        else:

            if horizon_point2[1] > horizon_point1[1]: #horizon part of a DOES NOT cross critical point

                surf_Ba = np.argmin((AdS_distance((vecA1_AdS, vecA2_AdS)) + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS)) + AdS_distance(CP((horizon_point1, horizon_point2), rH)), AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance(CP((vecB1_AdS, vecB2_AdS), rH)) + AdS_distance((horizon_point1, horizon_point2)), AdS_distance(CP((vecA1_AdS, vecA2_AdS), rH)) + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS)) + AdS_distance((horizon_point1, horizon_point2)), AdS_distance((vecB1_AdS, vecB2_AdS)) + AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance(CP((horizon_point1, horizon_point2), rH))))

                if surf_Ba == 0:

                    Ba = v0 + v5 + v7 + v13

                elif surf_Ba == 1:

                    Ba = v4 + v6 + v3 + v12

                elif surf_Ba ==2:

                    Ba = v1 + v5 + v7 + v12

                else:

                    Ba = v2 + v4 + v6 + v13

            else: # horizon_point2[1] < horizon_point1[1]: #horizon part of a crosses critical point

                surf_Ba = np.argmin((AdS_distance((vecA1_AdS, vecA2_AdS)) + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS)) + AdS_distance((horizon_point1, horizon_point2)), AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance(CP((vecB1_AdS, vecB2_AdS), rH)) + AdS_distance(CP((horizon_point1, horizon_point2), rH)), AdS_distance(CP((vecA1_AdS, vecA2_AdS), rH)) + AdS_distance((geod_point2, vecB1_AdS)) + AdS_distance((geod_point1, vecB2_AdS)) + AdS_distance(CP((horizon_point1, horizon_point2), rH)), AdS_distance((vecB1_AdS, vecB2_AdS)) + AdS_distance((vecA1_AdS, geod_point1)) + AdS_distance((vecA2_AdS, geod_point2)) + AdS_distance((horizon_point1, horizon_point2))))

                if surf_Ba == 0:

                    Ba = v0 + v5 + v7 + v12

                elif surf_Ba == 1:

                    Ba = v4 + v6 + v3 + v13

                elif surf_Ba == 2:

                    Ba = v1 + v5 + v7 + v13

                else:

                    Ba = v2 + v4 + v6 + v12

            #ABa------------------------------------------------------------

        if horizon_point1 == horizon_point2:

            if horizon_point1[1] < 0: #horizon belongs to a

                ABa = v5 + v7

            else: # horizon_point1[1] >= 0: #horizon belongs to b

                ABa = v5 + v7 + v12 + v13

        else:

            if horizon_point2[1] > horizon_point1[1]: #horizon part of b crosses critical point

                ABa = v5 + v7 + v13

            else: # horizon_point2[1] < horizon_point1[1]: #horizon part of b DOES NOT cross critical point

                ABa = v5 + v7 + v12


        Bab = A

        Aab = B

        ABb = a

        ab = AB

        Bb = Aa

        Ab = Ba

        b = ABa



        #I_AB-----------------------------------

        #surf_CM = A + B - AB 

        #E_sq-----------------------------------

        #surf_CM = Aa + Ba - ABa - a 

        #E_P------------------------------------

        #surf_CM = Aa 

        #E_R------------------------------------

        #surf_CM = AB + 2*Aa - a - ABa

        #E_Q------------------------------------

        surf_CM = A + B + Aa - Ba

        #surf_CM = 

        return surf_CM


    #Draw the surfaces

    surf_CM = surf_CM(E_CM.x)
    #surf_CM = [0,0,1,0,1,1,1,1,0,0,0,0,0,1]
    #print(surf_CM)

    color = (0, 'green', 'red')
    thickness = (0,2,6,9,12,15)

    if surf_CM[0] != 0:

        ax.plot(draw_BTZ_geod((vecA1_AdS, vecA2_AdS), rH)[0],draw_BTZ_geod((vecA1_AdS, vecA2_AdS), rH)[1], color = color[np.sign(surf_CM[0])], linewidth = 2)

    if surf_CM[1] != 0:

        ax.plot(draw_BTZ_geod(CP((vecA1_AdS, vecA2_AdS), rH), rH)[0],draw_BTZ_geod(CP((vecA1_AdS, vecA2_AdS), rH), rH)[1], color = color[np.sign(surf_CM[1])], linewidth = 2)

    if surf_CM[2] != 0:

        ax.plot(draw_BTZ_geod((vecB1_AdS, vecB2_AdS), rH)[0],draw_BTZ_geod((vecB1_AdS, vecB2_AdS), rH)[1], color = color[np.sign(surf_CM[2])], linewidth = 2)

    if surf_CM[3] != 0:

        ax.plot(draw_BTZ_geod(CP((vecB1_AdS, vecB2_AdS), rH), rH)[0],draw_BTZ_geod(CP((vecB1_AdS, vecB2_AdS), rH), rH)[1], color = color[np.sign(surf_CM[3])], linewidth = 2)

    if surf_CM[4] != 0:

        ax.plot(draw_BTZ_geod((vecA2_AdS, min_points_C[1]), rH)[0],draw_BTZ_geod((vecA2_AdS, min_points_C[1]), rH)[1], color = color[np.sign(surf_CM[4])], linewidth = 2)#, linestyle = '--')

    if surf_CM[5] != 0:

        ax.plot(draw_BTZ_geod((vecB1_AdS, min_points_C[1]), rH)[0],draw_BTZ_geod((vecB1_AdS, min_points_C[1]), rH)[1], color = color[np.sign(surf_CM[5])], linewidth = 2)

    if surf_CM[6] != 0:

        ax.plot(draw_BTZ_geod((vecA1_AdS, min_points_C[0]), rH)[0],draw_BTZ_geod((vecA1_AdS, min_points_C[0]), rH)[1], color = color[np.sign(surf_CM[6])], linewidth = 2)#, linestyle = '--')

    if surf_CM[7] != 0:

        ax.plot(draw_BTZ_geod((vecB2_AdS, min_points_C[0]), rH)[0],draw_BTZ_geod((vecB2_AdS, min_points_C[0]), rH)[1], color = color[np.sign(surf_CM[7])], linewidth = 2)

    if surf_CM[8] != 0:

        ax.plot(draw_BTZ_geod((min_points_C[0], min_points_C[1]), rH)[0],draw_BTZ_geod((min_points_C[0], min_points_C[1]), rH)[1], color = color[np.sign(surf_CM[8])], linewidth = 2)

    if surf_CM[9] != 0:

        ax.plot(draw_BTZ_geod(CP((min_points_C[0], min_points_C[1]), rH), rH)[0],draw_BTZ_geod(CP((min_points_C[0], min_points_C[1]), rH), rH)[1], color = color[np.sign(surf_CM[9])], linewidth = 2)

    if surf_CM[10] != 0:

        ax.plot(draw_BTZ_geod((min_points_C[0], min_points_C[2]), rH)[0],draw_BTZ_geod((min_points_C[0], min_points_C[2]), rH)[1], color = color[np.sign(surf_CM[10])], linewidth = 2)

    if surf_CM[11] != 0:

        ax.plot(draw_BTZ_geod((min_points_C[1], min_points_C[3]), rH)[0],draw_BTZ_geod((min_points_C[1], min_points_C[3]), rH)[1], color = color[np.sign(surf_CM[11])], linewidth = 2)

    if min_points_C[2] == min_points_C[3]:

        if surf_CM[12] == surf_CM[13] and surf_CM[12] != 0:

            ax.plot(draw_BTZ_geod(((0, -hor_max), (0, hor_max)), rH)[0],draw_BTZ_geod(((0, -hor_max), (0, hor_max)), rH)[1], color = color[np.sign(surf_CM[12])], linewidth = 2)

        elif surf_CM[12] != surf_CM[13]:

            print('horizon error?')

    else:

        if surf_CM[12] != 0:

            ax.plot(draw_BTZ_geod((min_points_C[2], min_points_C[3]), rH)[0],draw_BTZ_geod((min_points_C[2], min_points_C[3]), rH)[1], color = color[np.sign(surf_CM[12])], linewidth = 2)

        if surf_CM[13] != 0:

            ax.plot(draw_BTZ_geod(CP((min_points_C[2], min_points_C[3]), rH), rH)[0],draw_BTZ_geod(CP((min_points_C[2], min_points_C[3]), rH), rH)[1], color = color[np.sign(surf_CM[13])], linewidth = 2)

    

    #ax.plot(draw_BTZ_geod(CP((vecA1_AdS, vecB2_AdS), rH), rH)[0],draw_BTZ_geod(CP((vecA1_AdS, vecB2_AdS), rH), rH)[1], color = 'green', linewidth = 2)

    
    #FINISH DRAWING------------------------------------

    #--------------------------------------------------
    #-----THIS PART FOR FIGURES------------------------
    #--------------------------------------------------      

    surf_CM2 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    #print(surf_CM2)

    color = (0, 'green', 'red')
    thickness = (0,1,6,9,12,15)

    if surf_CM2[0] != 0:

        ax.plot(draw_BTZ_geod((vecA1_AdS, vecA2_AdS), rH)[0],draw_BTZ_geod((vecA1_AdS, vecA2_AdS), rH)[1], color = color[np.sign(surf_CM2[0])], linewidth = 2, linestyle = '--')

    if surf_CM2[1] != 0:

        ax.plot(draw_BTZ_geod(CP((vecA1_AdS, vecA2_AdS), rH), rH)[0],draw_BTZ_geod(CP((vecA1_AdS, vecA2_AdS), rH), rH)[1], color = color[np.sign(surf_CM2[1])], linewidth = 2, linestyle = '--')

    if surf_CM2[2] != 0:

        ax.plot(draw_BTZ_geod((vecB1_AdS, vecB2_AdS), rH)[0],draw_BTZ_geod((vecB1_AdS, vecB2_AdS), rH)[1], color = color[np.sign(surf_CM2[2])], linewidth = 2, linestyle = '--')

    if surf_CM2[3] != 0:

        ax.plot(draw_BTZ_geod(CP((vecB1_AdS, vecB2_AdS), rH), rH)[0],draw_BTZ_geod(CP((vecB1_AdS, vecB2_AdS), rH), rH)[1], color = color[np.sign(surf_CM2[3])], linewidth = 2, linestyle = '--')

    if surf_CM2[4] != 0:

        ax.plot(draw_BTZ_geod((vecA2_AdS, min_points_C[1]), rH)[0],draw_BTZ_geod((vecA2_AdS, min_points_C[1]), rH)[1], color = color[np.sign(surf_CM2[4])], linewidth = 2, linestyle = '--')

    if surf_CM2[5] != 0:

        ax.plot(draw_BTZ_geod((vecB1_AdS, min_points_C[1]), rH)[0],draw_BTZ_geod((vecB1_AdS, min_points_C[1]), rH)[1], color = color[np.sign(surf_CM2[5])], linewidth = 5, linestyle = (0,(15,15)))

    if surf_CM2[6] != 0:

        ax.plot(draw_BTZ_geod((vecA1_AdS, min_points_C[0]), rH)[0],draw_BTZ_geod((vecA1_AdS, min_points_C[0]), rH)[1], color = color[np.sign(surf_CM2[6])], linewidth = 2, linestyle = '--')

    if surf_CM2[7] != 0:

        ax.plot(draw_BTZ_geod((vecB2_AdS, min_points_C[0]), rH)[0],draw_BTZ_geod((vecB2_AdS, min_points_C[0]), rH)[1], color = color[np.sign(surf_CM2[7])], linewidth = 5, linestyle = (0,(15,15)))

    if surf_CM2[8] != 0:

        ax.plot(draw_BTZ_geod((min_points_C[0], min_points_C[1]), rH)[0],draw_BTZ_geod((min_points_C[0], min_points_C[1]), rH)[1], color = color[np.sign(surf_CM2[8])], linewidth = 2, linestyle = '--')

    if surf_CM2[9] != 0:

        ax.plot(draw_BTZ_geod(CP((min_points_C[0], min_points_C[1]), rH), rH)[0],draw_BTZ_geod(CP((min_points_C[0], min_points_C[1]), rH), rH)[1], color = color[np.sign(surf_CM2[9])], linewidth = 2, linestyle = '--')

    if surf_CM2[10] != 0:

        ax.plot(draw_BTZ_geod((min_points_C[0], min_points_C[2]), rH)[0],draw_BTZ_geod((min_points_C[0], min_points_C[2]), rH)[1], color = color[np.sign(surf_CM2[10])], linewidth = 2, linestyle = '--')

    if surf_CM2[11] != 0:

        ax.plot(draw_BTZ_geod((min_points_C[1], min_points_C[3]), rH)[0],draw_BTZ_geod((min_points_C[1], min_points_C[3]), rH)[1], color = color[np.sign(surf_CM2[11])], linewidth = 2, linestyle = '--')

    if min_points_C[2] == min_points_C[3]:

        if surf_CM2[12] == surf_CM2[13] and surf_CM2[12] != 0:

            ax.plot(draw_BTZ_geod(((0, -hor_max), (0, hor_max)), rH)[0],draw_BTZ_geod(((0, -hor_max), (0, hor_max)), rH)[1], color = color[np.sign(surf_CM2[12])], linewidth = 2, linestyle = '--')

        elif surf_CM2[12] != surf_CM2[13]:

            print('horizon error?')

    else:

        if surf_CM2[12] != 0:

            ax.plot(draw_BTZ_geod((min_points_C[2], min_points_C[3]), rH)[0],draw_BTZ_geod((min_points_C[2], min_points_C[3]), rH)[1], color = color[np.sign(surf_CM2[12])], linewidth = 2, linestyle = '--')

        if surf_CM2[13] != 0:

            ax.plot(draw_BTZ_geod(CP((min_points_C[2], min_points_C[3]), rH), rH)[0],draw_BTZ_geod(CP((min_points_C[2], min_points_C[3]), rH), rH)[1], color = color[np.sign(surf_CM2[13])], linewidth = 3, linestyle = (0,(15,15)))

    #linestyle = (0,(3,3))

    arc = pat.Arc((0,0), 1.8*rH, 1.8*rH, angle = 0, theta1 = 90, theta2 = 270, linewidth = 2, color = 'green')
    ax.add_artist(arc)
    

    #--------------------------------------------
    #--------------------------------------------
    #--------------------------------------------

    ax.plot(x_A1,y_A1,'o', color = 'black')
    ax.plot(x_B1,y_B1,'o', color = 'black')
    ax.plot(x_A2,y_A2,'o', color = 'black')
    ax.plot(x_B2,y_B2,'o', color = 'black')
    ax.plot(gp1_BTZ[0],gp1_BTZ[1],'o', color = 'black')
    ax.plot(gp2_BTZ[0],gp2_BTZ[1],'o', color = 'black')
    ax.plot(hp1_BTZ[0],hp1_BTZ[1],'o', color = 'black')
    ax.plot(hp2_BTZ[0],hp2_BTZ[1],'o', color = 'black')
    

    plt.savefig('path.filename')#, transparent = True)
            
#Bi_Corr_Meas(-8*m.pi/12, 8*m.pi/12, 9*m.pi/12, -9*m.pi/12, 0.3)
Bi_Corr_Meas(-15*m.pi/40, 15*m.pi/40, 25*m.pi/40, -25*m.pi/40, 0.27)
#Bi_Corr_Meas(-2.36, 2.36, 2.62, -2.62, 0.1)

#for rH in (0.1, 0.3, 0.5):
    
    #for B_size in range(1,11):

        #for A_size in range(1,12):

            #for A_center in range(0,12):

                #theta_A1 = (A_center - A_size)*m.pi/12
                #theta_A2 = (A_center + A_size)*m.pi/12
                #theta_B1 = (12 - B_size)*m.pi/12
                #theta_B2 = (-12 + B_size)*m.pi/12

                #if theta_A2 < theta_B1:

				    #print('In units of pi/12, A_size = ', A_size, ',   B_size = ', B_size, ',   A_center = ', A_center)
				    #print(str(round(theta_A1,2)) + '_' + str(round(theta_A2,2)) + '_' + str(round(theta_B1,2)) + '_' + str(round(theta_B2,2)) + '_' + str(rH))
                    #Bi_Corr_Meas(theta_A1, theta_A2, theta_B1, theta_B2, rH)
                


                    

            
            
            












                
        



