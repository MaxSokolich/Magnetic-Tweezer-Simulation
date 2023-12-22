#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 17:36:07 2021

@author: bizzaro
"""


#MATCH GAUSS CONTOURS IN PAPER TO theSE--->

#FIGURE OUT CONTOUR LINES -------->



import numpy as np
import matplotlib.pyplot as plt

gauss = (10**(-4)) #Tesla

milli = (10**(-3))
micro = (10**(-6))
nano = (10**(-9))
pico = (10**(-12))
#MAgnetic Monopole Approximation

'''
a point magneti charge associated iwth the magnetic 
flux can be defined similar to an electric charge as
'''

Pole_Pole_Distance = 405 * 2 * micro
Pole_Tip_Radius = 30 * micro

perm0 = ((4*np.pi)*(10**(-7)))

km = perm0/(4*np.pi)

reluctance = (1.8 * (10**9)) #A/Wb matched from ANSYS simulation
'''
#cross section area for reluctance from pole to pole through air is hard to calcualte since the magentic field spreads out through the air, 
so this value was determine by matching the reluctance found from ANSYS
'''
A = 0.18 *milli*milli #np.pi*(Pole_Tip_Radius**2) -----------> 0.18 mm^2 was found by matching ansys simulation
analystic_reluctance = (Pole_Pole_Distance/2)/(perm0 * (A))
print("Calculated Reluctance = ", analystic_reluctance *nano, "x 10^9")

#magnetomotive force = N*I
magnetomotive_force_left = 9
magnetomotive_force_top = 0
magnetomotive_force_right = 3
magnetomotive_force_bottom = 0


#Measuremtn Location
xx = -200*micro
yy = 0*micro
'''
phi = magnetomotive_force/reluctance #magnetc flux/ calculted thru magnetic curcuit analysis --- Weber = T/m^2 = kg/s^2/A  - ---  B * A

q = phi/perm0 #A/m^2

magnetomotive_force = phi * reluctance
'''



 

#R1= np.array([405*micro, 0*micro, 0*micro]) #unit vector point from location of magnetic charge to location in workspace
#r1 = np.linalg.norm(R1)
#u1 =  R1/r1

#B1 = (km * q1/(r1**2)) * u1


#For Q1

#R1 = np.array([-405*micro, -405*micro]) - np.array([x, y])
#r1 = np.linalg.norm(R1)
#u1 =  R1/r1


'''
Analytical Magnetic Monopole Approximation
'''



Left_Pole = [-Pole_Pole_Distance/2,0]
phi_left = magnetomotive_force_left/reluctance
Q_Left = phi_left/perm0

Top_Pole = [0,Pole_Pole_Distance/2]
phi_top = magnetomotive_force_top/reluctance
Q_Top = phi_top/perm0

Right_Pole = [Pole_Pole_Distance/2,0]
phi_right = magnetomotive_force_right/reluctance
Q_Right = phi_right/perm0

Bottom_Pole = [0,-Pole_Pole_Distance/2]
phi_bottom = magnetomotive_force_bottom/reluctance
Q_Bottom = phi_bottom/perm0



'''
the magnetic field produce by the magnetic charge in teh workspace is 
'''   

#Bottom Pole
Bx1 = ((km * (Q_Bottom/(np.sqrt((xx-(0))**2 +(yy-(-(Pole_Pole_Distance/2)))**2)**2))) * ((xx-0)/np.sqrt((xx-(0))**2 +(yy-(-(Pole_Pole_Distance/2)))**2)))
By1 = ((km * (Q_Bottom/(np.sqrt((xx-(0))**2 +(yy-(-(Pole_Pole_Distance/2)))**2)**2))) * ((yy-(-(Pole_Pole_Distance/2)))/np.sqrt((xx-(0))**2 +(yy-(-(Pole_Pole_Distance/2)))**2)))

#Top Pole
Bx2 = ((km * (Q_Top/(np.sqrt((xx-(0))**2 +(yy-((Pole_Pole_Distance/2)))**2)**2))) * ((xx-0)/np.sqrt((xx-(0))**2 +(yy-((Pole_Pole_Distance/2)))**2)))
By2 = ((km * (Q_Top/(np.sqrt((xx-(0))**2 +(yy-((Pole_Pole_Distance/2)))**2)**2))) * ((yy-((Pole_Pole_Distance/2)))/np.sqrt((xx-(0))**2 +(yy-((Pole_Pole_Distance/2)))**2)))

#Left Pole
Bx3 = ((km * (Q_Left/(np.sqrt((xx-(-(Pole_Pole_Distance/2)))**2 +(yy-(0))**2)**2))) * ((xx-(-(Pole_Pole_Distance/2)))/np.sqrt((xx-(-(Pole_Pole_Distance/2)))**2 +(yy-(0))**2)))
By3 = ((km * (Q_Left/(np.sqrt((xx-(-(Pole_Pole_Distance/2)))**2 +(yy-(0))**2)**2))) * ((yy-(0))/np.sqrt((xx-(-(Pole_Pole_Distance/2)))**2 +(yy-(0))**2)))

#Right Pole
Bx4 = ((km * (Q_Right/(np.sqrt((xx-((Pole_Pole_Distance/2)))**2 +(yy-(0))**2)**2))) * ((xx-((Pole_Pole_Distance/2)))/np.sqrt((xx-((Pole_Pole_Distance/2)))**2 +(yy-(0))**2)))
By4 = ((km * (Q_Right/(np.sqrt((xx-((Pole_Pole_Distance/2)))**2 +(yy-(0))**2)**2))) * ((yy-(0))/np.sqrt((xx-((Pole_Pole_Distance/2)))**2 +(yy-(0))**2)))

#Thru super position
Bx = Bx1 + Bx2 + Bx3  + Bx4
By = By1 + By2 + By3  + By4

BTot = [1000*Bx,1000*By]
print("BVector = ", BTot)



#STREAMPLOT CONTROUR PLOT
X = np.linspace(-(Pole_Pole_Distance/2), (Pole_Pole_Distance/2), 4)
Y = np.linspace(-(Pole_Pole_Distance/2),(Pole_Pole_Distance/2), 4)
x,y = np.meshgrid(X, Y)



def Bfield_X(q,Pole_Location):
    x1 = Pole_Location[0]
    y1 = Pole_Location[1]
    
    r = (np.sqrt((x-x1)**2 + (y-y1)**2))
    
    ux = ((x-x1)/r)
    
    Bx = km * (q/(r**2)) * ux
    
    return Bx


def Bfield_Y(q,Pole_Location):
    x1 = Pole_Location[0]
    y1 = Pole_Location[1]
    
    r = (np.sqrt((x-x1)**2 + (y-y1)**2))
      
    uy = ((y-y1)/r)
    
    By = km * (q/(r**2)) * uy
    
    return By 


#Thru Super Positon

BTotal_X = (1/gauss)*(Bfield_X(Q_Left,Left_Pole) + Bfield_X(Q_Top,Top_Pole) + Bfield_X(Q_Right,Right_Pole) +  Bfield_X(Q_Bottom,Bottom_Pole))
BTotal_Y = (1/gauss)*(Bfield_Y(Q_Left,Left_Pole) + Bfield_Y(Q_Top,Top_Pole) + Bfield_Y(Q_Right,Right_Pole) +  Bfield_Y(Q_Bottom,Bottom_Pole))

BVect_Gauss = [BTotal_X ,BTotal_Y]


plt.figure()

plt.xlim(-(Pole_Pole_Distance/2),(Pole_Pole_Distance/2))
plt.ylim(-(Pole_Pole_Distance/2),(Pole_Pole_Distance/2))
plt.xlabel("X (m)")
plt.ylabel("Y (m)")

speed = np.sqrt((BTotal_X)**2+(BTotal_Y)**2)

plt.streamplot(x,y,BTotal_X,BTotal_Y, color=speed)
plt.figure()
cs = plt.contourf(x,y,speed,cmap='coolwarm')

lines = []
for line in cs.collections[3].get_paths():
    lines.append(line.vertices)
    
#plt.plot(lines[0][:, 0], lines[0][:, 1],linewidth = 5,color="k")
plt.colorbar()
plt.show()
