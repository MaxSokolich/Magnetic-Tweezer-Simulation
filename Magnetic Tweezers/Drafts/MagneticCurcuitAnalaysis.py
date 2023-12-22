#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 09:11:14 2021

@author: bizzaro
"""

'''
4 Pole Magnetic Curcuit Analysis

F = magnetic motive force = N*I
phi = magnetic flux = B*Area
reluctance = L/(perm0*permr*Area)
'''
import numpy as np

milli = (10**(-3))
micro = (10**(-6))
nano = (10**(-9))
pico = (10**(-12))

perm0 = ((4*np.pi)*(10**(-7)))
Pole_Pole_Distance = 405 * 2 * micro

'''
In this circuit, Ra is significantly larger than Rp and Ry 
as the yoke and poles are made of materials much more magnetic permeable than air.
'''
#Yoke Reluctance
yoke_length = (25 *milli)
YpermR = 5000
Yoke_AreaC = (25*milli*milli)
Ry = (yoke_length)/(YpermR*perm0*Yoke_AreaC)
print("Yoke Reluctance = ", Ry)
#Pole Reluctance

pole_length = (17 *milli)
PpermR = 3600
Pole_AreaC = (0.5*milli*milli)
Rp = (pole_length)/(PpermR*perm0*Pole_AreaC)
print("Pole Reluctance = ", Rp)

#Air Gap Reluctance (Length measured from pole to workspace center)
Air_length = (Pole_Pole_Distance/2)
ApermR = 1
Air_AreaC = (0.18*milli*milli) 
'''
#-------IMPORTANT--------This is what should be calculted for as              
the relutance of the air gap was calcualted from ANSYS simulation. It is hard to calculate this because
the magnetic field is spread out 

Thus the model will fail if pole diameter is 2 big for 2 reasons;  monopole approximation doesnt hold, 
it will be harder to model the magnetic curcuit and obtain the relationship between input current 
and magnetic flux
'''
Ra = (Air_length)/(ApermR*perm0*Air_AreaC)
print("Air Gap Reluctance = ", Ra)


#HOPKINSONS LAW 
'''
basically equivlant to electrical crucit ohms law but for magnetic curcuits
magnetomotive force = phi * reluctance
voltage =  current * resistance

The magnetic flux through each pole can
therefore be determined according to the magnetic circuit. where pole and yoke reluctances
are neglected due to 4 orders of magnetiudde less signifcant than air gap
'''

#phi_matrix = np.array([phi1,phi2,ph3,phi4])
I1 = 0
I2 = 5
I3 = 0
I4 = 0

NTurns = 21
KI = np.array([[3/4, -1/4, -1/4, -1/4],  #matrix KI is the distribution matrix of the magnetic flux.
              [-1/4, 3/4, -1/4, -1/4],
              [-1/4, -1/4, 3/4, -1/4],
              [-1/4, -1/4, -1/4, 3/4]])

I_matrix = np.array([I1,I2,I3,I4])

Q_Matric_Const = (NTurns/(perm0*Ra))
Qmatrix = Q_Matric_Const * (KI @ I_matrix)

phi1 = Qmatrix[0]
phi2 = Qmatrix[1]
phi3 = Qmatrix[2]
phi4 = Qmatrix[3]
print("magnetic flux from pole1 = ", phi1, "Wb")
print("magnetic flux from pole2 = ", phi2, "Wb")
print("magnetic flux from pole3 = ", phi3, "Wb")
print("magnetic flux from pole4 = ", phi4, "Wb")

#Q is a column matrix



