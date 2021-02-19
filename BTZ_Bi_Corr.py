#!/usr/bin/env python3
import numpy as np
import math as m
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as pat
import scipy.optimize as opt
import random as rand
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

def PoincareGeodesicDistance(vec1, vec2):
	Deltax = vec2[0]-vec1[0]
	Deltay = vec2[1]-vec1[1]
	r1 = m.sqrt(vec1[0]**2 + vec1[1]**2)
	r2 = m.sqrt(vec2[0]**2 + vec2[1]**2)
	dist = m.acosh(1 + 2*(Deltax**2 + Deltay**2)/((1 - r1**2)*(1-r2**2)))
	#print(vec1,vec2)
	#print(dist)
	return(m.acosh(1 + 2*(Deltax**2 + Deltay**2)/((1 - r1**2)*(1-r2**2))))



#define a function which takes two polar points in BTZ and gives the angular coordinate phi_1 which achieves the SHORTEST geodesic connecting them

#CartesianBTZDiskPoints = [PolarToCartesian(BTZDiskPoints[0]), PolarToCartesian(BTZDiskPoints[1])] 
#PoincareDiskPoints = [] 
#for vec in BTZDiskPoints: 

 #   PoincareDiskPoints.append(BTZDiskToPoincareDisk(vec, BTZDiskToBTZ(rH)))    
#CartesianPoincareDiskPoints = [] 
#for vec in PoincareDiskPoints: 
#    CartesianPoincareDiskPoints.append(PolarToCartesian(vec)) 

#def BTZ_min_coords(vec1, vec2, rH):
	
	##The following is BTZ Disk points in polar coordinates
	#BTZ_polar = (vec1, vec2)
	#BTZ_polar1 = ((vec1[0],vec1[1] + 2*m.pi), vec2)
	#BTZ_polar2 = ((vec1[0],vec1[1] - 2*m.pi), vec2)

	##print(BTZ_polar)
	##print(BTZ_polar1)
	##print(BTZ_polar2)

	#PoincareDiskPoints = (BTZDiskToPoincareDisk(BTZ_polar[0], BTZDiskToBTZ(rH)), BTZDiskToPoincareDisk(BTZ_polar[1], BTZDiskToBTZ(rH)))
	#PoincareDiskPoints1 = (BTZDiskToPoincareDisk(BTZ_polar1[0], BTZDiskToBTZ(rH)), BTZDiskToPoincareDisk(BTZ_polar1[1], BTZDiskToBTZ(rH)))
	#PoincareDiskPoints2 = (BTZDiskToPoincareDisk(BTZ_polar2[0], BTZDiskToBTZ(rH)), BTZDiskToPoincareDisk(BTZ_polar2[1], BTZDiskToBTZ(rH)))

	#CartesianPoincareDiskPoints =(PolarToCartesian(PoincareDiskPoints[0]),PolarToCartesian(PoincareDiskPoints[1]))
	#CartesianPoincareDiskPoints1 =(PolarToCartesian(PoincareDiskPoints1[0]),PolarToCartesian(PoincareDiskPoints1[1]))
	#CartesianPoincareDiskPoints2 =(PolarToCartesian(PoincareDiskPoints2[0]),PolarToCartesian(PoincareDiskPoints2[1]))

	#AdS_dist = PoincareGeodesicDistance(CartesianPoincareDiskPoints[0], CartesianPoincareDiskPoints[1])
	#AdS_dist1 = PoincareGeodesicDistance(CartesianPoincareDiskPoints1[0], CartesianPoincareDiskPoints1[1])
	#AdS_dist2 = PoincareGeodesicDistance(CartesianPoincareDiskPoints2[0], CartesianPoincareDiskPoints2[1])
	
	##print(AdS_dist,AdS_dist1,AdS_dist2)

	##AdS_polar = (BTZToPoincareDisk(BTZ_polar[0],BTZDiskToBTZ(rH)), BTZToPoincareDisk(BTZ_polar[1],BTZDiskToBTZ(rH)))
	##AdS_polar1 = (BTZToPoincareDisk(BTZ_polar1[0],BTZDiskToBTZ(rH)), BTZToPoincareDisk(BTZ_polar1[1],BTZDiskToBTZ(rH)))
	##AdS_polar2 = (BTZToPoincareDisk(BTZ_polar2[0],BTZDiskToBTZ(rH)), BTZToPoincareDisk(BTZ_polar2[1],BTZDiskToBTZ(rH)))

	##tryvec = (0.99999999, AdS_polar[0][1])
	##tryvec1 = (0.99999999, AdS_polar[1][1])

	##AdS_cart = (PolarToCartesian(tryvec),PolarToCartesian(tryvec1))
	##AdS_cart1 =a (PolarToCartesian(AdS_polar1[0]),PolarToCartesian(AdS_polar1[1]))
	##AdS_cart2 = (PolarToCartesian(AdS_polar2[0]),PolarToCartesian(AdS_polar2[1]))

	##AdS_dist = PoincareGeodesicDistance(AdS_cart[0], AdS_cart[1])
	##AdS_dist1 = PoincareGeodesicDistance(AdS_cart1[0], AdS_cart1[1])
	##AdS_dist2 = PoincareGeodesicDistance(AdS_cart2[0], AdS_cart2[1])

	#BTZ_min_geod = np.argmin((AdS_dist, AdS_dist1, AdS_dist2))

	##print(AdS_cart[0], AdS_cart[1])

	##print(BTZ_min_geod)

	#if BTZ_min_geod == 0:
		#nvec = BTZ_polar

	#elif BTZ_min_geod == 1:
		#nvec = BTZ_polar1

	#elif BTZ_min_geod == 2:
		#nvec = BTZ_polar2

	##print(nvec)
	##print('-----------------------------------------------------------------------')

	#return(nvec[0], nvec[1])

#BTZ_min_coords((0.9999999, 0.8),(0.99999999, 2.3),0.2)


#define a distance function which takes two polar points in BTZ and gives the length of the SHORTEST geodesic connecting them

def BTZ_distance(vec1, vec2, rH):

	coords = (vec1, vec2)	

	AdS_coords_polar = (BTZDiskToPoincareDisk(coords[0], BTZDiskToBTZ(rH)), BTZDiskToPoincareDisk(coords[1], BTZDiskToBTZ(rH)))

	AdS_cart = (PolarToCartesian(AdS_coords_polar[0]),PolarToCartesian(AdS_coords_polar[1]))

	BTZ_dist = PoincareGeodesicDistance(AdS_cart[0], AdS_cart[1])
	return BTZ_dist

Nmax = 50

def draw_BTZ_geod(vec1, vec2, rH):

	PoinGeodesicPoints = []

	if abs(vec1[0] - rH) < 0.0000000000001 and abs(vec2[0] - rH) < 0.0000000000001: #both points on horizon

		BTZDiskPoints = [(rH, vec1[1]), (rH, vec2[1])]

		PoincareDiskPoints = []

		for vec in BTZDiskPoints:
                        
			PoincareDiskPoints.append(BTZDiskToPoincareDisk(vec, BTZDiskToBTZ(rH)))

		CartesianPoinDiskPoints = []

		for vec in PoincareDiskPoints:
                        
			CartesianPoinDiskPoints.append(PolarToCartesian(vec))

		for n in range(0,Nmax):

			point = (0, min(CartesianPoinDiskPoints[0][1],CartesianPoinDiskPoints[1][1]) + (n/Nmax)*abs(CartesianPoinDiskPoints[0][1] - CartesianPoinDiskPoints[1][1]))

			PoinGeodesicPoints.append(point)      


	else: #at least one point not on horizon

		if abs(vec1[0] - rH) < 0.0000000000001: 

			BTZDiskPoints = ((rH, vec1[1]),vec2)

		elif abs(vec2[0] - rH) < 0.0000000000001:    

			BTZDiskPoints = (vec1, (rH, vec2[1]))

		else:

			BTZDiskPoints = (vec1, vec2)

		CartesianBTZDiskPoints = [PolarToCartesian(BTZDiskPoints[0]), PolarToCartesian(BTZDiskPoints[1])] 

		PoincareDiskPoints = [] #Red and Blue points in polar coordinates

		for vec in BTZDiskPoints: 

			PoincareDiskPoints.append(BTZDiskToPoincareDisk(vec, BTZDiskToBTZ(rH)))

		angle_1 = min(PoincareDiskPoints[0][1], PoincareDiskPoints[1][1])
		angle_2 = max(PoincareDiskPoints[0][1], PoincareDiskPoints[1][1])
		angle = angle_2 - angle_1

		CartesianPoincareDiskPoints = [] #Red and Blue points in carterisan coordinates

		for vec in PoincareDiskPoints: 

			CartesianPoincareDiskPoints.append(PolarToCartesian(vec))

                

         

                #Distance = PoincareGeodesicDistance(CartesianPoincareDiskPoints[0], CartesianPoincareDiskPoints[1]) 

		CircleData = PoincareGeodesicCircle(CartesianPoincareDiskPoints[0], CartesianPoincareDiskPoints[1]) 

		x_1 = CartesianPoincareDiskPoints[0][0]
		x_2 = CartesianPoincareDiskPoints[1][0]
		y_1 = CartesianPoincareDiskPoints[0][1]
		y_2 = CartesianPoincareDiskPoints[1][1]

		theta_ = m.atan2(y_2,x_2) - m.atan2(y_1,x_1)

		if theta_ > m.pi:

			X_1 = x_2
			Y_1 = y_2
			X_2 = x_1
			Y_2 = y_1

		elif theta_ < -m.pi:

			X_1 = x_1
			Y_1 = y_1
			X_2 = x_2
			Y_2 = y_2

		elif 0 < theta_ < m.pi:
                        
			X_1 = x_1
			Y_1 = y_1
			X_2 = x_2
			Y_2 = y_2
		else:

			X_1 = x_2
			Y_1 = y_2
			X_2 = x_1
			Y_2 = y_1
                
		c_x = CircleData[0][0]
		c_y = CircleData[0][1]

		theta_2 = m.atan2(Y_2 - c_y, X_2 - c_x)
		theta_1 = m.atan2(Y_1 - c_y, X_1 - c_x)

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
   

#Black line is the boundary, red line is the horizon. 

#Green lines are phi = plus/minus pi on the BTZ disk and are identified. 

#Area contained within red and green lines and segment of the boundary is one fundamental region of the BTZ disk. 

#Purple line is the geodesic connecting the two points. 

#It may be useful to have one of the two points outside the fundamental region to find the shortest geodesic. 


#start with two non-overlapping boundary regions and an optimized correlation measure

cutoff = 0.0000001
max_r = 1 - cutoff

def Bi_Corr_Meas(theta_A1, theta_A2, theta_B1, theta_B2, rH):

	vecA1 = (max_r, theta_A1)
	vecA2 = (max_r, theta_A2)
	vecB1 = (max_r, theta_B1)
	vecB2 = (max_r, theta_B2)	

	x_A1 = m.cos(theta_A1)*max_r
	y_A1 = m.sin(theta_A1)*max_r
	x_A2 = m.cos(theta_A2)*max_r
	y_A2 = m.sin(theta_A2)*max_r

	x_B1 = m.cos(theta_B1)*max_r
	y_B1 = m.sin(theta_B1)*max_r
	x_B2 = m.cos(theta_B2)*max_r
	y_B2 = m.sin(theta_B2)*max_r


	#make a figure

	fig, axs = plt.subplots(2,3)

	for CM in range(1,6):

		if CM == 1:

			SP = [0,1]

		if CM == 2:

			SP = [0,2]

		if CM == 3:

			SP = [1,0]

		if CM == 4:

			SP = [1,1]

		if CM == 5:

			SP = [0,0]

		axs[SP[0],SP[1]].set_xlim((-1,1))
		axs[SP[0],SP[1]].set_ylim((-1,1))
	
		axs[SP[0],SP[1]].set_aspect('equal')
		
		axs[SP[0],SP[1]].tick_params(top = 'off', bottom = 'off', left = 'off', right = 'off', labeltop = 'off', labelbottom = 'off', labelleft = 'off', labelright = 'off' )

		#draw boundary, boundary region endpoints, and horizon

		boundary = plt.Circle((0,0), 1, color = 'black', fill = False, clip_on = False)
		axs[SP[0],SP[1]].add_artist(boundary)

		arc1 = pat.Arc((0,0), 2, 2, angle = 0, theta1 = 180*theta_A1/m.pi, theta2 = 180*theta_A2/m.pi, linewidth = 3, color = 'blue')
		arc2 = pat.Arc((0,0), 2, 2, angle = 0, theta1 = 180*theta_B1/m.pi, theta2 = 180*theta_B2/m.pi, linewidth = 3, color = 'blue')
		axs[SP[0],SP[1]].add_artist(arc1)
		axs[SP[0],SP[1]].add_artist(arc2)


		#horizon = plt.Circle((0,0), rH, color = 'black', fill = False, clip_on = False)
		#axs[SP[0],SP[1]].add_artist(horizon)	

		#HorizonPointsPolar = []
		#HorizonPointsCart = []

		#Nmax = 500

		#for n in range(1, Nmax):

			#x = rH*m.cos(2*m.pi*n/Nmax)

			#y = rH*m.sin(2*m.pi*n/Nmax)

			#HorizonPointsPolar.append((rH, 2*m.pi*n/Nmax))
			#HorizonPointsCart.append((x,y))

		#axs[SP[0],SP[1]].scatter(*zip(*HorizonPointsCart),s=5, color = 'purple') 



		#plot boundary region endpoints

		axs[SP[0],SP[1]].plot((x_A1), (y_A1), 'o', color = 'black')
		axs[SP[0],SP[1]].plot((x_A2), (y_A2), 'o', color = 'black')	
		axs[SP[0],SP[1]].plot((x_B1), (y_B1), 'o', color = 'black')
		axs[SP[0],SP[1]].plot((x_B2), (y_B2), 'o', color = 'black')

		c_1 = BTZ_distance(vecA2, vecB1, rH) + BTZ_distance(vecB2, vecA1, rH) + BTZ_distance((rH,m.pi),(rH,-m.pi),rH)
		c_2 = BTZ_distance(vecA1, vecA2, rH) + BTZ_distance(vecB1, (max_r, theta_B2 + 2*m.pi), rH)
		c_3 = BTZ_distance(vecA2, vecB1, rH) + BTZ_distance((max_r, theta_B2 + 2*m.pi), vecA1, rH)
		c_4 = BTZ_distance(vecB2, vecA1, rH) + BTZ_distance(vecB1, (max_r, theta_A2 + 2*m.pi), rH)

		print(c_1,c_2,c_3,c_4)

		EW_type = np.argmin((c_1, c_2, c_3, c_4))

		S_A = min(BTZ_distance(vecA1, vecA2, rH), BTZ_distance((max_r, theta_A1 + 2*m.pi), vecA2, rH) + BTZ_distance((rH, -m.pi), (rH, m.pi), rH))
		S_B = min(BTZ_distance(vecB1, (max_r, theta_B2 + 2*m.pi), rH), BTZ_distance(vecB1, vecB2, rH) + BTZ_distance((rH, -m.pi), (rH, m.pi), rH))
		S_AB = min(c_1, c_2, c_3, c_4)

		I_AB = S_A + S_B - S_AB

		#print('I(A:B) = ', I_AB)

		#print(EW_options)
		#print(EW_type)

		#draw entanglement wedge

		if EW_type == 0:
		
			axs[SP[0],SP[1]].scatter(*zip(*draw_BTZ_geod(vecA1, vecB2, rH)), s=1, color = 'purple') 
			axs[SP[0],SP[1]].scatter(*zip(*draw_BTZ_geod(vecA2, vecB1, rH)), s=1, color = 'purple')
			axs[SP[0],SP[1]].scatter(*zip(*draw_BTZ_geod((rH, -m.pi), (rH, m.pi), rH)), s=1, color = 'purple')

		elif EW_type == 1:

			print('Entanglement wedge disconnected')
			return()

		elif EW_type == 2 or EW_type == 3:

			print('Entanglement wedge connected but does not contain black hole')
			return()

		#Define objective function to be minimized

		
		
		def f_CM(opt_points):

			geod_point1 = CartesianToPolar(draw_BTZ_geod(vecA1, vecB2, rH)[opt_points[0]])
			geod_point2 = CartesianToPolar(draw_BTZ_geod(vecA2, vecB1, rH)[opt_points[1]])

			horizon_point1 = CartesianToPolar(draw_BTZ_geod((rH, -m.pi), (rH, m.pi), rH)[opt_points[2]])
			horizon_point2 = CartesianToPolar(draw_BTZ_geod((rH, -m.pi), (rH, m.pi), rH)[opt_points[3]])

			#S_A--------------------------------------------------------------

			#defined above

			#S_B--------------------------------------------------------------

			#defined above

			#S_a--------------------------------------------------------------

			if horizon_point1 == horizon_point2:

				if horizon_point1[1] < 0: #horizon belongs to a

					S_a = BTZ_distance(vecA1, geod_point1, rH) + BTZ_distance(vecA2, geod_point2, rH) + BTZ_distance((rH, -m.pi), (rH, m.pi), rH)

				elif horizon_point1[1] >= 0: #horizon belongs to b

					S_a = BTZ_distance(vecA1, geod_point1, rH) + BTZ_distance(vecA2, geod_point2, rH)

			else:

				if horizon_point2[1] > horizon_point1[1]: #horizon part of a DOES NOT cross critical point

					S_a = BTZ_distance(vecA1, geod_point1, rH) + BTZ_distance(vecA2, geod_point2, rH) + BTZ_distance(horizon_point1, horizon_point2, rH)

				elif horizon_point2[1] < horizon_point1[1]: #horizon part of a crosses critical point

					S_a = BTZ_distance(vecA1, geod_point1, rH) + BTZ_distance(vecA2, geod_point2, rH) + BTZ_distance((rH, horizon_point2[1] + 2*m.pi), horizon_point1, rH)

			#S_AB-------------------------------------------------------------

			#defined above

			#S_Aa-------------------------------------------------------------

			if horizon_point1 == horizon_point2:

				if horizon_point1[1] < 0: #horizon belongs to a

					S_Aa = min(BTZ_distance(geod_point1, geod_point2, rH) + BTZ_distance((rH, -m.pi), (rH, m.pi), rH), BTZ_distance((geod_point1[0], geod_point1[1] + 2*m.pi), geod_point2, rH))

				elif horizon_point1[1] >= 0: #horizon belongs to b

					
					S_Aa = min(BTZ_distance(geod_point1, geod_point2, rH), BTZ_distance((geod_point1[0], geod_point1[1] + 2*m.pi), geod_point2, rH) + BTZ_distance((rH, -m.pi), (rH, m.pi), rH))

			else:

				S_Aa = BTZ_distance(geod_point1, horizon_point1, rH) + BTZ_distance(geod_point2, horizon_point2, rH)

			#S_Ba-------------------------------------------------------------

			if horizon_point1 == horizon_point2:

				if horizon_point1[1] < 0: #horizon belongs to a

					S_Ba = min(BTZ_distance(vecA1, vecA2, rH) + BTZ_distance(geod_point2, vecB1, rH) + BTZ_distance(geod_point1, vecB2, rH), BTZ_distance((max_r, theta_A1 + 2*m.pi), vecA2, rH) + BTZ_distance((rH, -m.pi), (rH, m.pi), rH) + BTZ_distance(geod_point2, vecB1, rH) + BTZ_distance(geod_point1, vecB2, rH), BTZ_distance(vecA1, geod_point1, rH) + BTZ_distance(vecA2, geod_point2, rH) + BTZ_distance(vecB1, (max_r, theta_B2 + 2*m.pi), rH) + BTZ_distance((rH, -m.pi), (rH, m.pi), rH), BTZ_distance(vecA1, geod_point1, rH) + BTZ_distance(vecA2, geod_point2, rH) + BTZ_distance(vecB1, vecB2, rH))

				elif horizon_point1[1] >= 0: #horizon belongs to b

					S_Ba = min(BTZ_distance(vecA1, vecA2, rH) + BTZ_distance(geod_point2, vecB1, rH) + BTZ_distance(geod_point1, vecB2, rH) + BTZ_distance((rH, -m.pi), (rH, m.pi), rH), BTZ_distance((max_r, theta_A1 + 2*m.pi), vecA2, rH) + BTZ_distance(geod_point2, vecB1, rH) + BTZ_distance(geod_point1, vecB2, rH), BTZ_distance(vecA1, geod_point1, rH) + BTZ_distance(vecA2, geod_point2, rH) + BTZ_distance(vecB1, (max_r, theta_B2 + 2*m.pi), rH), BTZ_distance(vecA1, geod_point1, rH) + BTZ_distance(vecA2, geod_point2, rH) + BTZ_distance(vecB1, vecB2, rH) + BTZ_distance((rH, -m.pi), (rH, m.pi), rH))

			else:

				if horizon_point2[1] > horizon_point1[1]: #horizon part of a DOES NOT cross critical point

					S_Ba = min(BTZ_distance(vecA1, vecA2, rH) + BTZ_distance(geod_point2, vecB1, rH) + BTZ_distance(geod_point1, vecB2, rH) + BTZ_distance(horizon_point2, (rH, horizon_point1[1] + 2*m.pi), rH), BTZ_distance(vecA1, geod_point1, rH) + BTZ_distance(vecA2, geod_point2, rH) + BTZ_distance(vecB1, (max_r, theta_B2 + 2*m.pi), rH) + BTZ_distance(horizon_point1, horizon_point2, rH))

				elif horizon_point2[1] < horizon_point1[1]: #horizon part of a crosses critical point

					S_Ba = min(BTZ_distance(vecA1, vecA2, rH) + BTZ_distance(geod_point2, vecB1, rH) + BTZ_distance(geod_point1, vecB2, rH) + BTZ_distance(horizon_point1, horizon_point2, rH), BTZ_distance(vecA1, geod_point1, rH) + BTZ_distance(vecA2, geod_point2, rH) + BTZ_distance(vecB1, (max_r, theta_B2 + 2*m.pi), rH) + BTZ_distance(horizon_point1, (rH, horizon_point2[1] + 2*m.pi), rH))

			#S_ABa------------------------------------------------------------

			if horizon_point1 == horizon_point2:

				if horizon_point1[1] < 0: #horizon belongs to a

					S_ABa = BTZ_distance(vecB1, geod_point2, rH) + BTZ_distance(vecB2, geod_point1, rH)

				elif horizon_point1[1] >= 0: #horizon belongs to b

					S_ABa = BTZ_distance(vecB1, geod_point2, rH) + BTZ_distance(vecB2, geod_point1, rH) + BTZ_distance((rH, -m.pi), (rH, m.pi), rH)

			else:

				if horizon_point2[1] > horizon_point1[1]: #horizon part of b crosses critical point

					S_ABa = BTZ_distance(vecB1, geod_point2, rH) + BTZ_distance(vecB2, geod_point1, rH) + BTZ_distance((rH, horizon_point1[1] + 2*m.pi), horizon_point2, rH)

				elif horizon_point2[1] < horizon_point1[1]: #horizon part of b DOES NOT cross critical point

					S_ABa = BTZ_distance(vecB1, geod_point2, rH) + BTZ_distance(vecB2, geod_point1, rH) + BTZ_distance(horizon_point1, horizon_point2, rH)

		
			S_Bab = S_A

			S_Aab = S_B

			S_ABb = S_a

			S_ab = S_AB

			S_Bb = S_Aa

			S_Ab = S_Ba

			S_b = S_ABa

		

			if CM == 1: #E_P

				f_CM = S_Aa
			
			elif CM == 2: #E_Q

				f_CM = S_A + S_B + S_Aa - S_Ba

			elif CM == 3: #E_sq

				f_CM = S_Aa + S_Ba - S_ABa - S_a

			elif CM == 4: #E_R

				f_CM = S_AB + 2*S_Aa - S_a - S_ABa
			
			elif CM == 5: #I_AB

				f_CM = S_A + S_B - S_AB

			return f_CM

		#Minimize the objective

		E_CM = 10000

		#low resolution search

		if CM != 5:

			#print(len(draw_BTZ_geod(vecA1, vecB2, rH)),len(draw_BTZ_geod(vecA2, vecB1, rH)),len(draw_BTZ_geod((rH, -m.pi), (rH, m.pi), rH)),len(draw_BTZ_geod((rH, -m.pi), (rH, m.pi), rH)))
			#A = len(draw_BTZ_geod(vecA1, vecB2, rH))*len(draw_BTZ_geod(vecA2, vecB1, rH))*len(draw_BTZ_geod((rH, -m.pi), (rH, m.pi), rH))*len(draw_BTZ_geod((rH, -m.pi), (rH, m.pi), rH))
			#print(A/10000)

			for geod_point1 in range(0, len(draw_BTZ_geod(vecA1, vecB2, rH)),10):

				for geod_point2 in range(0,len(draw_BTZ_geod(vecA2, vecB1, rH)),10):

					for horizon_point1 in range(0,len(draw_BTZ_geod((rH, -m.pi), (rH, m.pi), rH)),10):

						for horizon_point2 in range(0,len(draw_BTZ_geod((rH, -m.pi), (rH, m.pi), rH)),10):
							print(geod_point1,geod_point2,horizon_point1,horizon_point2)

							if f_CM((geod_point1, geod_point2, horizon_point1, horizon_point2)) < E_CM:

								E_CM = f_CM((geod_point1, geod_point2, horizon_point1, horizon_point2))
								min_points_low = [geod_point1, geod_point2, horizon_point1, horizon_point2]

						
	

			#For the endpoints: 	


			if min_points_low[0] + 6 > len(draw_BTZ_geod(vecA1, vecB2, rH)) - 1:

				min_points_low[0] = min_points_low[0] - len(draw_BTZ_geod(vecA1, vecB2, rH))

			if min_points_low[1] + 6 > len(draw_BTZ_geod(vecA2, vecB1, rH)) - 1:

				min_points_low[1] = min_points_low[1] - len(draw_BTZ_geod(vecA2, vecB1, rH))


			if min_points_low[2] + 6 > len(draw_BTZ_geod((rH, -m.pi), (rH, m.pi), rH)) - 1:

				min_points_low[2] = min_points_low[2] - len(draw_BTZ_geod((rH, -m.pi), (rH, m.pi), rH))

			if min_points_low[3] + 6 > len(draw_BTZ_geod((rH, -m.pi), (rH, m.pi), rH)) - 1:

				min_points_low[3] = min_points_low[3] - len(draw_BTZ_geod((rH, -m.pi), (rH, m.pi), rH))
	

			#high resolution search

			for geod_point1 in range(min_points_low[0] - 5, min_points_low[0] + 6):

				for geod_point2 in range(min_points_low[1] - 5, min_points_low[1] + 6):

					for horizon_point1 in range(min_points_low[2] -5, min_points_low[2] + 6):

						for horizon_point2 in range(min_points_low[3] - 5, min_points_low[3] + 6):
							print(geod_point1,geod_point2,horizon_point1,horizon_point2)

							if f_CM((geod_point1, geod_point2, horizon_point1, horizon_point2)) <= E_CM:

								E_CM = f_CM((geod_point1, geod_point2, horizon_point1, horizon_point2))
								min_points_ = (geod_point1, geod_point2, horizon_point1, horizon_point2)

			#print('E_CM = ', E_CM)
			#print('I(A:B) = ', I_AB)
			#print('S(A) = ', S_A)
			#print('S(B) = ', S_B)
			#print(min_points_)

			geod_point_1 = CartesianToPolar(draw_BTZ_geod(vecA1, vecB2, rH)[min_points_[0]])
			geod_point_2 = CartesianToPolar(draw_BTZ_geod(vecA2, vecB1, rH)[min_points_[1]])
			horizon_point_1 = CartesianToPolar(draw_BTZ_geod((rH, -m.pi), (rH, m.pi), rH)[min_points_[2]])
			horizon_point_2 = CartesianToPolar(draw_BTZ_geod((rH, -m.pi), (rH, m.pi), rH)[min_points_[3]])

			min_points = (geod_point_1, geod_point_2, horizon_point_1, horizon_point_2)
			#print(min_points)

			#Define a function which returns a vector giving the bulk surfaces contributing to f_CM (order = (A1-A2, A1-A2 CP, B1-B2, B1-B2 CP, A2-G2, B1-G2, A1-G1, B2-G1, G1-G2, G1-G2 CP, G1-H1, G2-H2, H1-H2, H1-H2 CP)) 

		def surf_CM(x):

			geod_point1 = CartesianToPolar(draw_BTZ_geod(vecA1, vecB2, rH)[x[0]])
			geod_point2 = CartesianToPolar(draw_BTZ_geod(vecA2, vecB1, rH)[x[1]])
			horizon_point1 = CartesianToPolar(draw_BTZ_geod((rH, -m.pi), (rH, m.pi), rH)[x[2]])
			horizon_point2 = CartesianToPolar(draw_BTZ_geod((rH, -m.pi), (rH, m.pi), rH)[x[3]])

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

			surf_A = np.argmin((BTZ_distance(vecA1, vecA2, rH), BTZ_distance((max_r, theta_A1 + 2*m.pi), vecA2, rH) + BTZ_distance((rH, -m.pi), (rH, m.pi), rH)))

			if surf_A == 0:

				A = v0

			else:

				A = v1 + v12 + v13

			#B--------------------------------------------------------------

			surf_B = np.argmin((BTZ_distance(vecB1, (max_r, theta_B2 + 2*m.pi), rH), BTZ_distance(vecB1, vecB2, rH) + BTZ_distance((rH, -m.pi), (rH, m.pi), rH)))

			if surf_B == 0:

				B = v3

			else:

				B = v2 + v12 + v13

			#a--------------------------------------------------------------

			if horizon_point1 == horizon_point2:

				if horizon_point1[1] < 0: #horizon belongs to a

					a = v4 + v6 + v12 + v13

				elif horizon_point1[1] >= 0: #horizon belongs to b

					a = v4 + v6

			else:

				if horizon_point2[1] > horizon_point1[1]: #horizon part of a DOES NOT cross critical point

					a = v4 + v6 + v12

				elif horizon_point2[1] < horizon_point1[1]: #horizon part of a crosses critical point

					a = v4 + v6 + v13

			#AB-------------------------------------------------------------

			AB = v6 + v7 + v4 + v5 + v12 + v13

			#Aa-------------------------------------------------------------

			if horizon_point1 == horizon_point2:

				if horizon_point1[1] < 0: #horizon belongs to a

					surf_Aa = np.argmin((BTZ_distance(geod_point1, geod_point2, rH) + BTZ_distance((rH, -m.pi), (rH, m.pi), rH), BTZ_distance((geod_point1[0], geod_point1[1] + 2*m.pi), geod_point2, rH)))

					if surf_Aa == 0:

						Aa = v8 + v12 + v13

					else:

						Aa = v9

				elif horizon_point1[1] >= 0: #horizon belongs to b

					
					surf_Aa = np.argmin((BTZ_distance(geod_point1, geod_point2, rH), BTZ_distance((geod_point1[0], geod_point1[1] + 2*m.pi), geod_point2, rH) + BTZ_distance((rH, -m.pi), (rH, m.pi), rH)))

					if surf_Aa == 0:

						Aa = v8

					else:

						Aa = v9 + v12 + v13

			else:

				Aa = v10 + v11

			#Ba-------------------------------------------------------------

			if horizon_point1 == horizon_point2:

				if horizon_point1[1] < 0: #horizon belongs to a

					surf_Ba = np.argmin((BTZ_distance(vecA1, vecA2, rH) + BTZ_distance(geod_point2, vecB1, rH) + BTZ_distance(geod_point1, vecB2, rH), BTZ_distance((max_r, theta_A1 + 2*m.pi), vecA2, rH) + BTZ_distance((rH, -m.pi), (rH, m.pi), rH) + BTZ_distance(geod_point2, vecB1, rH) + BTZ_distance(geod_point1, vecB2, rH), BTZ_distance(vecA1, geod_point1, rH) + BTZ_distance(vecA2, geod_point2, rH) + BTZ_distance(vecB1, (max_r, theta_B2 + 2*m.pi), rH) + BTZ_distance((rH, -m.pi), (rH, m.pi), rH), BTZ_distance(vecA1, geod_point1, rH) + BTZ_distance(vecA2, geod_point2, rH) + BTZ_distance(vecB1, vecB2, rH)))

					if surf_Ba == 0:

						Ba = v0 + v5 + v7

					elif surf_Ba == 1:

						Ba = v1 + v12 + v13 + v5 + v7

					elif surf_Ba == 2:

						Ba = v4 + v6 + v3 + v12 + v13

					else:

						Ba = v4 + v6 + v2

				elif horizon_point1[1] >= 0: #horizon belongs to b

					surf_Ba = np.argmin((BTZ_distance(vecA1, vecA2, rH) + BTZ_distance(geod_point2, vecB1, rH) + BTZ_distance(geod_point1, vecB2, rH) + BTZ_distance((rH, -m.pi), (rH, m.pi), rH), BTZ_distance((max_r, theta_A1 + 2*m.pi), vecA2, rH) + BTZ_distance(geod_point2, vecB1, rH) + BTZ_distance(geod_point1, vecB2, rH), BTZ_distance(vecA1, geod_point1, rH) + BTZ_distance(vecA2, geod_point2, rH) + BTZ_distance(vecB1, (max_r, theta_B2 + 2*m.pi), rH), BTZ_distance(vecA1, geod_point1, rH) + BTZ_distance(vecA2, geod_point2, rH) + BTZ_distance(vecB1, vecB2, rH) + BTZ_distance((rH, -m.pi), (rH, m.pi), rH)))

					if surf_Ba == 0:

						Ba = v0 + v5 + v7 + v12 + v13

					elif surf_Ba == 1:

						Ba = v1 + v5 + v7

					elif surf_Ba == 2:

						Ba = v4 + v6 + v3

					else:

						Ba = v4 + v6 + v2 + v12 + v13

			else:

				if horizon_point2[1] > horizon_point1[1]: #horizon part of a DOES NOT cross critical point

					surf_Ba = np.argmin((BTZ_distance(vecA1, vecA2, rH) + BTZ_distance(geod_point2, vecB1, rH) + BTZ_distance(geod_point1, vecB2, rH) + BTZ_distance(horizon_point2, (rH, horizon_point1[1] + 2*m.pi), rH), BTZ_distance(vecA1, geod_point1, rH) + BTZ_distance(vecA2, geod_point2, rH) + BTZ_distance(vecB1, (max_r, theta_B2 + 2*m.pi), rH) + BTZ_distance(horizon_point1, horizon_point2, rH)))

					if surf_Ba == 0:

						Ba = v0 + v5 + v7 + v13

					else:

						Ba = v4 + v6 + v3 + v12

				elif horizon_point2[1] < horizon_point1[1]: #horizon part of a crosses critical point

					surf_Ba = np.argmin((BTZ_distance(vecA1, vecA2, rH) + BTZ_distance(geod_point2, vecB1, rH) + BTZ_distance(geod_point1, vecB2, rH) + BTZ_distance(horizon_point1, horizon_point2, rH), BTZ_distance(vecA1, geod_point1, rH) + BTZ_distance(vecA2, geod_point2, rH) + BTZ_distance(vecB1, (max_r, theta_B2 + 2*m.pi), rH) + BTZ_distance(horizon_point1, (rH, horizon_point2[1] + 2*m.pi), rH)))

					if surf_Ba == 0:

						Ba = v0 + v5 + v7 + v12

					else:

						Ba = v4 + v6 + v3 + v13

			#ABa------------------------------------------------------------

			if horizon_point1 == horizon_point2:

				if horizon_point1[1] < 0: #horizon belongs to a

					ABa = v5 + v7

				elif horizon_point1[1] >= 0: #horizon belongs to b

					ABa = v5 + v7 + v12 + v13

			else:

				if horizon_point2[1] > horizon_point1[1]: #horizon part of b crosses critical point

					ABa = v5 + v7 + v13

				elif horizon_point2[1] < horizon_point1[1]: #horizon part of b DOES NOT cross critical point

					ABa = v5 + v7 + v12

		
			Bab = A

			Aab = B

			ABb = a

			ab = AB

			Bb = Aa

			Ab = Ba

			b = ABa

		

			if CM == 1: #E_P

				surf_CM = Aa

			elif CM == 2: #E_Q

				surf_CM = A + B + Aa - Ba

			elif CM == 3: #E_sq

				surf_CM = Aa + Ba - ABa - a

			elif CM == 4: #E_R

				surf_CM = AB + 2*Aa - a - ABa

			elif CM == 5: #I_AB

				surf_CM = A + B - AB

			#surf_CM = A + B + Bb - Ab

			return surf_CM


		#Draw the surfaces

		if CM != 5:

			surf_CM = surf_CM(min_points_)
			print(surf_CM)

		elif CM == 5:

			surf_CM = surf_CM((0,0,0,0))

		color = (0, 'green', 'red')
		thickness = (0,3,6,9,12,15)

		if surf_CM[0] != 0:

			axs[SP[0],SP[1]].scatter(*zip(*draw_BTZ_geod(vecA1, vecA2, rH)), s = thickness[abs(surf_CM[0])], color = color[np.sign(surf_CM[0])])

		if surf_CM[1] != 0:

			axs[SP[0],SP[1]].scatter(*zip(*draw_BTZ_geod((max_r, theta_A1 + 2*m.pi), vecA2, rH)), s = thickness[abs(surf_CM[1])], color = color[np.sign(surf_CM[1])])

		if surf_CM[2] != 0:

			axs[SP[0],SP[1]].scatter(*zip(*draw_BTZ_geod(vecB1, vecB2, rH)), s = thickness[abs(surf_CM[2])], color = color[np.sign(surf_CM[2])])	

		if surf_CM[3] != 0:

			axs[SP[0],SP[1]].scatter(*zip(*draw_BTZ_geod(vecB1, (max_r, theta_B2 + 2*m.pi), rH)), s = thickness[abs(surf_CM[3])], color = color[np.sign(surf_CM[3])])

		if surf_CM[4] != 0:

			axs[SP[0],SP[1]].scatter(*zip(*draw_BTZ_geod(vecA2, min_points[1], rH)), s = thickness[abs(surf_CM[4])], color = color[np.sign(surf_CM[4])])

		if surf_CM[5] != 0:

			axs[SP[0],SP[1]].scatter(*zip(*draw_BTZ_geod(vecB1, min_points[1], rH)), s = thickness[abs(surf_CM[5])], color = color[np.sign(surf_CM[5])])

		if surf_CM[6] != 0:

			axs[SP[0],SP[1]].scatter(*zip(*draw_BTZ_geod(vecA1, min_points[0], rH)), s = thickness[abs(surf_CM[6])], color = color[np.sign(surf_CM[6])])

		if surf_CM[7] != 0:

			axs[SP[0],SP[1]].scatter(*zip(*draw_BTZ_geod(vecB2, min_points[0], rH)), s = thickness[abs(surf_CM[7])], color = color[np.sign(surf_CM[7])])

		if surf_CM[8] != 0:

			axs[SP[0],SP[1]].scatter(*zip(*draw_BTZ_geod(min_points[0], min_points[1], rH)), s = thickness[abs(surf_CM[8])], color = color[np.sign(surf_CM[8])])

		if surf_CM[9] != 0:

			axs[SP[0],SP[1]].scatter(*zip(*draw_BTZ_geod((min_points[0][0], min_points[0][1] + 2*m.pi), min_points[1], rH)), s = thickness[abs(surf_CM[9])], color = color[np.sign(surf_CM[9])])

		if surf_CM[10] != 0:

			axs[SP[0],SP[1]].scatter(*zip(*draw_BTZ_geod(min_points[0], min_points[2], rH)), s = thickness[abs(surf_CM[10])], color = color[np.sign(surf_CM[10])])

		if surf_CM[11] != 0:

			axs[SP[0],SP[1]].scatter(*zip(*draw_BTZ_geod(min_points[1], min_points[3], rH)), s = thickness[abs(surf_CM[11])], color = color[np.sign(surf_CM[11])])

		if min_points_[2] == min_points_[3]:

			if surf_CM[12] == surf_CM[13] and surf_CM[12] != 0:

				axs[SP[0],SP[1]].scatter(*zip(*draw_BTZ_geod((rH, -m.pi), (rH, m.pi), rH)), s = thickness[abs(surf_CM[12])], color = color[np.sign(surf_CM[12])])

			elif surf_CM[12] != surf_CM[13]:

				print('horizon error?')

		else:

			if surf_CM[12] != 0:

				axs[SP[0],SP[1]].scatter(*zip(*draw_BTZ_geod(min_points[2], min_points[3], rH)), s = thickness[abs(surf_CM[12])], color = color[np.sign(surf_CM[12])])

			if surf_CM[13] != 0:

				axs[SP[0],SP[1]].scatter(*zip(*draw_BTZ_geod((rH, min(min_points[2][1], min_points[3][1]) + 2*m.pi), (rH, max(min_points[2][1], min_points[3][1])), rH)), s = thickness[abs(surf_CM[13])], color = color[np.sign(surf_CM[13])])

		#FINISH DRAWING--------------------------------------------------

		axs[1,2].axis('off')

		if CM == 1:

			name = 'P'
			E = str(round(E_CM, 4))

		if CM == 2:

			name = 'Q'
			E = str(round(E_CM, 4))

		if CM == 3:

			name = 'sq'
			E = str(round(E_CM, 4))

		if CM == 4:

			name = 'R'
			E = str(round(E_CM, 4))

		if CM == 5:

			name = 'I'
			E = str(round(I_AB, 4))

		S_minpoints = '''({G1_r}, {G1_phi}), ({G2_r}, {G2_phi})
phi1 = {H1_phi}, phi2 = {H2_phi}'''.format(G1_r = str(round(min_points[0][0],2)), G1_phi = str(round(min_points[0][1],2)), G2_r = str(round(min_points[1][0],2)), G2_phi = str(round(min_points[1][1],2)), H1_phi = str(round(min_points[2][1],2)), H2_phi = str(round(min_points[3][1],2)))

		S = '''E_{CM} = {E_CM}
{surf}
{min_points}'''.format(CM = name, E_CM = E, surf = str(surf_CM), min_points = S_minpoints)

		axs[SP[0],SP[1]].set_xlabel(S)

		#if min_points[2] = min_points[3]:

		#	if CM == 1:

		#		axs[SP[0],SP[1]].scatter(*zip(*draw_BTZ_geod(min_points[0], min_points[1], rH)), s=5, color = 'red')

		#else:

		#	if CM == 1:

		#		axs[SP[0],SP[1]].scatter(*zip(*draw_BTZ_geod(min_points[0], min_points[2], rH)), s=5, color = 'red')
		#		axs[SP[0],SP[1]].scatter(*zip(*draw_BTZ_geod(min_points[1], min_points[3], rH)), s=5, color = 'red')

		#delta_h = BTZ_distance(min_points[2], min_points[3], rH) - BTZ_distance((rH, min(min_points[2][1], min_points[3][1]) + 2*m.pi), (rH, max(min_points[2][1], min_points[3][1])), rH)
		#delta_g = BTZ_distance(vecA1, vecA2, rH) - BTZ_distance(vecB1, (max_r, theta_B2 + 2*m.pi), rH)
		#print(delta_h, delta_g)	
		#plt.title(str(CM) + ', ' + str(theta_A1) + ', ' + str(theta_A2) + ', ' + str(theta_B1) + ', ' + str(theta_B2) + ', ' + str(rH) + ', E_CM = ' + str(E_CM))

				

		#------------END OF CM LOOP-------------------------------



	S = '''theta_A1 = {theta_A1}
theta_A2 = {theta_A2}
theta_B1 = {theta_B1}
theta_B2 = {theta_B2}
rH = {rH}
S(A) = {S_A}
S(B) = {S_B}'''.format(theta_A1 = str(round(theta_A1,2)), theta_A2 = str(round(theta_A2,2)), theta_B1 = str(round(theta_B1,2)), theta_B2 = str(round(theta_B2,2)), rH = str(rH), S_A = str(round(S_A,2)), S_B = str(round(S_B,2)))
	fig.text(0.9, 0.12, S, verticalalignment='bottom', horizontalalignment='right')
	fig.tight_layout()
	disp = (theta_A2 + (theta_B1 - theta_A2)/2 - m.pi/2)/(m.pi/24)
	S_disp = str(round(disp))
	plt.savefig('/home/jole1266/Figures/'+ S_disp + '_' + str(round(theta_A1,2)) + '_' + str(round(theta_A2,2)) + '_' + str(round(theta_B1,2)) + '_' + str(round(theta_B2,2)) + '_' + str(rH) + '.png')

	return()


#Bi_Corr_Meas(4, -1.6, 1.6, m.pi - 1.2, -m.pi + 1.2, 0.5) 
#Bi_Corr_Meas(1, -1.1, 1.1, m.pi - 1.1, -m.pi + 1.1, 0.2)
#Bi_Corr_Meas(1, -0.8, 0.8, 0.9, 2.4, 0.5)
#Bi_Corr_Meas(1, -0.4, 1, m.pi - 1, -m.pi + 1.3, 0.2)

disp_angle = 0

for rH in (0.1, 0.3, 0.5): 

	for region_angle in range(1, 12):

		theta_A1 = -m.pi/2 - (disp_angle)*m.pi/24 + (region_angle)*m.pi/24
		theta_A2 = m.pi/2 + (disp_angle)*m.pi/24 - (region_angle)*m.pi/24
		theta_B1 = m.pi/2 + (disp_angle)*m.pi/24 + (region_angle)*m.pi/24
		theta_B2 = -m.pi/2 - (disp_angle)*m.pi/24 - (region_angle)*m.pi/24

		if theta_B1 < m.pi:

			Bi_Corr_Meas(theta_A1, theta_A2, theta_B1, theta_B2, rH)
	




