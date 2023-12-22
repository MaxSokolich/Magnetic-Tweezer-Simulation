#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 09:11:14 2021

@author: bizzaro
"""
'''
Analytical Magnetic Monopole Approximation


https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=5280242

F = magnetic motive force = N*I
phi = magnetic flux = B*Area
reluctance = L/(perm0*permr*Area)
'''

'''
phi = magnetomotive_force/reluctance #magnetc flux/ calculted thru magnetic curcuit analysis --- Weber = T/m^2 = kg/s^2/A  - ---  B * A

q = phi/perm0 #A/m^2

magnetomotive_force = phi * reluctance
'''



import numpy as np
import matplotlib.pyplot as plt

milli = (10**(-3))
micro = (10**(-6))
nano = (10**(-9))
pico = (10**(-12))

perm0 = ((4*np.pi)*(10**(-7)))
Pole_Pole_Distance = 405 * 2 * micro


km = perm0/(4*np.pi)

#Current inputs
ILeft = .9
ITop = -.6
IRight = .9
IBottom = -.6

NTurns = 21 #number of turns, assumed same for each pole coil

Air_AreaC = (0.18*milli*milli)  # cross sectioal area of 

Meas_Loc = [0*micro, 200*micro] #measure the magnetic field at a postion in the workspace


print("Gradientmax = ", "alot\n")

#MAGNETIC CURCUIT

'''
In this circuit, Ra is significantly larger than Rp and Ry 
as the yoke and poles are made of materials much more magnetic permeable than air.
'''
#Yoke Reluctance
yoke_length = (25 *milli)
YpermR = 5000
Yoke_AreaC = (25*milli*milli)
Ry = (yoke_length)/(YpermR*perm0*Yoke_AreaC)
#print("Yoke Reluctance = ", Ry)



#Pole Reluctance
pole_length = (17 *milli)
PpermR = 3600
Pole_AreaC = (0.5*milli*milli)
Rp = (pole_length)/(PpermR*perm0*Pole_AreaC)
#print("Pole Reluctance = ", Rp)

#Air Gap Reluctance (Length measured from pole to workspace center)
Air_length = (Pole_Pole_Distance/2)
ApermR = 1

'''
#-------IMPORTANT--------This is what should be calculted for as              
the relutance of the air gap was calcualted from ANSYS simulation. It is hard to calculate this because
the magnetic field is spread out 

Thus the model will fail if pole diameter is 2 big for 2 reasons;  monopole approximation doesnt hold, 
it will be harder to model the magnetic curcuit and obtain the relationship between input current 
and magnetic flux
'''
Ra = (Air_length)/(ApermR*perm0*Air_AreaC)
print("Air Gap Reluctance = ", Ra, "\n")









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

KI = np.array([[3/4, -1/4, -1/4, -1/4],  #matrix KI is the distribution matrix of the magnetic flux.
              [-1/4, 3/4, -1/4, -1/4],
              [-1/4, -1/4, 3/4, -1/4],
              [-1/4, -1/4, -1/4, 3/4]])

I_matrix = np.array([ILeft,ITop,IRight,IBottom])

print("Currents as L,T,R,B = ", I_matrix, "\n")
Q_Matric_Const = (NTurns/(perm0*Ra))
Qmatrix = Q_Matric_Const * (KI @ I_matrix)

'''
 constraint that the net magnetic flux is always
zero, as determined by the Gauss’s law for magnetism. According to the definition of magnetic charges, the summation of the
magnetic charges is thus equal to zero
'''

print("Q Charge Matrix = ", Qmatrix, "\n")
Q_Left = Qmatrix[0]
Q_Top = Qmatrix[1]
Q_Right = Qmatrix[2]
Q_Bottom = Qmatrix[3]

#Q is a column matrix


#MAGNETIC FIELD



Left_Pole = [-Pole_Pole_Distance/2,0]

Top_Pole = [0,Pole_Pole_Distance/2]

Right_Pole = [Pole_Pole_Distance/2,0]

Bottom_Pole = [0,-Pole_Pole_Distance/2]





def Bfield_X(q,Pole_Location,Meas_Loc):
    x1 = Pole_Location[0]
    y1 = Pole_Location[1]
    x = Meas_Loc[0]
    y = Meas_Loc[1]
    r = (np.sqrt((x-x1)**2 + (y-y1)**2))
    
    ux = ((x-x1)/r)
    
    Bx = km * (q/(r**2)) * ux
    
    return Bx


def Bfield_Y(q,Pole_Location, Meas_Loc):
    x1 = Pole_Location[0]
    y1 = Pole_Location[1]
    x = Meas_Loc[0]
    y = Meas_Loc[1]
    r = (np.sqrt((x-x1)**2 + (y-y1)**2))

    uy = ((y-y1)/r)
    
    By = km * (q/(r**2)) * uy
    
    return By 

BTX = (1000)*(Bfield_X(Q_Left,Left_Pole, Meas_Loc) + Bfield_X(Q_Top,Top_Pole, Meas_Loc) + Bfield_X(Q_Right,Right_Pole, Meas_Loc) +  Bfield_X(Q_Bottom,Bottom_Pole, Meas_Loc))
BTY = (1000)*(Bfield_Y(Q_Left,Left_Pole, Meas_Loc) + Bfield_Y(Q_Top,Top_Pole, Meas_Loc) + Bfield_Y(Q_Right,Right_Pole, Meas_Loc) +  Bfield_Y(Q_Bottom,Bottom_Pole, Meas_Loc))
BT = [BTX,BTY]
print("\nLocation = ",Meas_Loc)
print("Magnetic Field =", BT, 'mT\n')



#STREAMPLOT CONTROUR PLOT
X = np.linspace(-(Pole_Pole_Distance/2), (Pole_Pole_Distance/2), 4)
Y = np.linspace(-(Pole_Pole_Distance/2),(Pole_Pole_Distance/2), 4)
x,y = np.meshgrid(X, Y)


#Thru Super Positon

BTotal_X = (1000)*(Bfield_X(Q_Left,Left_Pole,[x,y]) + Bfield_X(Q_Top,Top_Pole,[x,y]) + Bfield_X(Q_Right,Right_Pole,[x,y]) +  Bfield_X(Q_Bottom,Bottom_Pole,[x,y]))
BTotal_Y = (1000)*(Bfield_Y(Q_Left,Left_Pole,[x,y]) + Bfield_Y(Q_Top,Top_Pole,[x,y]) + Bfield_Y(Q_Right,Right_Pole,[x,y]) +  Bfield_Y(Q_Bottom,Bottom_Pole,[x,y]))

BVect_Gauss = [BTotal_X ,BTotal_Y]


plt.figure()

plt.xlim(-(Pole_Pole_Distance/2),(Pole_Pole_Distance/2))
plt.ylim(-(Pole_Pole_Distance/2),(Pole_Pole_Distance/2))
plt.xlabel("X (m)")
plt.ylabel("Y (m)")

speed = np.sqrt((BTotal_X)**2+(BTotal_Y)**2)

plt.streamplot(x,y,BTotal_X,BTotal_Y, color=speed,cmap='coolwarm')
plt.figure()
cs = plt.contourf(x,y,speed,cmap='coolwarm')

lines = []
for line in cs.collections[3].get_paths():
    lines.append(line.vertices)
  
#plt.plot(lines[0][:, 0], lines[0][:, 1],linewidth = 5,color="k")
plt.colorbar()

plt.show()


