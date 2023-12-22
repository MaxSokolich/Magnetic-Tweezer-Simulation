#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 17:44:47 2021

@author: bizzaro

States 6mT at center

depth calcualtions page 38

"The saturation magnetization (Ms) is the maximum magnetic moment per unit volume for a magnetic material."
microospheres 

23.96 emu g−1.

or 0.8 T?


Do we want the input to be velocity aka distance traveled per time?

Or do we want the input a certain magnetic strength?





.8/100 slop of super paramagnetic sphere = perm = .008

.2/200 = .001 sloper for ferromagnetic sphere

Source =    https://www.spherotech.com/2020%20Product%20Detail%20Pages/Spherotech%20Paramagnetic%20Particles.pdf




"""

import numpy as np

milli = (10**(-3))
micro = (10**(-6))
nano = (10**(-9))
pico = (10**(-12))


def divergence(f):
    """
    Computes the divergence of the vector field f, corresponding to dFx/dx + dFy/dy + ...
    :param f: List of ndarrays, where every item of the list is one dimension of the vector field
    :return: Single ndarray of the same shape as each of the items in f, which corresponds to a scalar field
    """
    num_dims = len(f)
    return np.ufunc.reduce(np.add, [np.gradient(f[i], axis=i) for i in range(num_dims)])




#INPUTS
# %%
Pole_Pole_Distance = (1* milli)
Pole_Diameter = (40*micro)

bead_diameter = (10*micro)                                                                                                    #Bead Diameter in Microns
perm0 = ((4*np.pi)*10**(-7))                                                                                                  #Magnetic Vacuum Permeability
perm_bead = 40 
                                                                                                                              #Magnetic Microsphere Permeability
Viscosity =  (8.90 * (10**(-4)))      
                                                                                                                              #Fluid Viscosity Pa·s (Water)
Velocity = np.array([10*micro,0*micro,0*micro])                                                                                           #Desired Microsphere Velocity [Vx,Vy,Vz] in meters/s

N_Turns = 527                                                                                                                 #Coil Number of turns      
'''
Magnetic charge inputs. distances form each actuated pole to the microbot
These would be calculated from camera system
'''
#Pole 1
Vector_1 = np.array([500*micro,0*micro,0*micro])                                                                                #vector direction pointing from pole 1 to the magnetic microrobot
U_Vector_1 = Vector_1/np.linalg.norm(Vector_1)  
R_1 = np.linalg.norm(Vector_1)                                                                                                #the distance from pole 1 tip to the magnetic microrobot

#Pole 2                                                                                                                       #the distance from pole 2 tip to the magnetic microrobot
Vector_2 = np.array([500*micro,0*micro,0*micro])                                                                                #vector direction pointing from pole 2 to the magnetic microrobot
U_Vector_2 = Vector_2/np.linalg.norm(Vector_2)  
R_2 = np.linalg.norm(Vector_2)                                                                                                #the distance from pole 1 tip to the magnetic microrobot

#Pole 3                                                                                                                       #the distance from pole 3 tip to the magnetic microrobot
Vector_3 = np.array([500*micro,0*micro,0*micro])                                                                                #vector direction pointing from pole 3 to the magnetic microrobot
U_Vector_3 = Vector_3/np.linalg.norm(Vector_3)  
R_3 = np.linalg.norm(Vector_3)                                                                                                #the distance from pole 1 tip to the magnetic microrobot


'''
Ignoring flux leakage, the reluctance of 
the circuit is simply the sum of the reluctances of the yoke (y), 
pole tips (p) and the air gap (a) between the poles.
Because the permeability and average cross-section of the yoke 
and pole tips are several orders of magnitude higher than that of the air gap,
the total reluctance of the magnetic circuit is determined only by the air gap. 
(Pg 64 Vries thesis)
'''

Reluctance = Pole_Pole_Distance/(perm0 * ((np.pi*Pole_Diameter**2)/4))                                                        #reluctance between the pole tip and the working space center


'''
#direction vector, Xm, in the measurement coordinate system given from the 
user or the closed loop control algorithm
'''
Measurment_Direction = np.array([40*micro,40*micro,0*micro]).transpose()                                                       #Direction Vector (distance and direction) to point microbead in (measured from measurement coord system or camera coordiante system)
# %%
#
#
#
#
#
#
#
#
#
#
#COORDINATE SYSTEM TRANSFORMATION
# %%

'''
The transformation matrix from the measurement coordinate system to the 
actuation coordinate system is computed.
Rx and Ry are the rotation matrices for space-fixed coordinate 
transformation around x-axis and y-axis, respectively. 
When the direction vector, "Measurment_Direction ", in the measurement coordinate 
system is given from the user or the closed loop control algorithm, 
"Transformation_Matrix" will transform it into vector "Acutation_Direction " in the 
actuation coordinate system
'''


Transformation_Matrix = [[0.8165, 0, -0.5774], 
                         [0.4082, 0.7071, 0.5774], 
                         [0.4082, -0.7071, 0.5774]]



Acutation_Direction = Measurment_Direction#np.dot(Transformation_Matrix, Measurment_Direction)                                #Direction Vector to point microbead along (measured from measurment actuation coordinate system)
Actuation_Unit_Direction = Acutation_Direction / np.linalg.norm(Acutation_Direction)                                          #Normalized Actuation direction
print("Actuation Unit Direction = ",Actuation_Unit_Direction, "\n")
# %%
#
#
#
#
#
#
#
#
#
#

#CURRENT CALCULATIONS
# %%

Amplifier = 1        
                                                                                                                              #Amplitude of gain parameter to output adequate current
Compensation_Vector = np.array([0.04,0.01,0.01])                                                                              #compensation vector for current input is determined through experimental calibration,compensates for unequal coils

Final_Current_Output = Amplifier * Actuation_Unit_Direction + Compensation_Vector                                             #The actual current amplitude sent to the electronics (the sign of each index will determine which opposite facing coil gets power)

print("Final Current Output = ",Final_Current_Output, "A\n")
# %%
#
#
#
#
#
#
#
#
#
#

# MAGNETIC FIELD CALCULATIONS
# %%
'''
The magnetic field generated from three poles can be expressed.
B is the magnetic flux density, or amount of magnetic force induced on the microbot due to a magnetizing force H.
ri is the distance from the magnetic pole tip to the magnetic microrobot, 
μ0 is the permeability of free space (vacuum), 
qi is the magnetic charge defined by q = Φ/μ0, 
Φ is the magnetic flux, 
ui is the unit vector pointing from the magnetic charge to the magnetic microrobot

'''

  
Magnetic_Charges = (N_Turns/(perm0*Reluctance)) * Final_Current_Output



Magnetic_Charge_1 = Magnetic_Charges[0]                                                                                       #magnetic charge generated from pole 1
Magnetic_Charge_2 = Magnetic_Charges[1]                                                                                       #magnetic charge generated from pole 2
Magnetic_Charge_3 = Magnetic_Charges[2]                                                                                       #magnetic charge generated from pole 3

Magnetic_Flux_1 = ((4*np.pi)/perm0)  *  ((Magnetic_Charge_1)/(R_1**2))  *  U_Vector_1                                         #Magnetic flux generated at the microbots location from pole 1
Magnetic_Flux_2 = ((4*np.pi)/perm0)  *  ((Magnetic_Charge_2)/(R_2**2))  *  U_Vector_2                                         #Magnetic field generated at the microbots location from pole 2
Magnetic_Flux_3 = ((4*np.pi)/perm0)  *  ((Magnetic_Charge_3)/(R_3**3))  *  U_Vector_3                                         #Magnetic field generated at the microbots location from pole 3

Magnetic_Flux_Density = Magnetic_Flux_1  +  Magnetic_Flux_2  +  Magnetic_Flux_3
print("Magnetic Flux Density From Pole 1 = ",Magnetic_Flux_1, "T")
print("Magnetic Flux Density From Pole 2 = ",Magnetic_Flux_2, "T")
print("Magnetic Flux Density From Pole 3 = ",Magnetic_Flux_3, "T\n")

print("Total Magnetic Flux Density on Microbot = ", Magnetic_Flux_Density,"T\n\n")
# %%
#
#
#
#
#
#
#
#
#
#

#MAGNETIC FORCE CALCUATIONS
# %%
"""
This is the gradient of the magnetic field, B [T], acting on a particle with a magnetic mo- ment, m [Am2]. 
The latter can be estimated from its volume integral, i.e. the bulk mass magnetization, M [Am2/kg or emu/g],
m=VM

https://arxiv.org/pdf/1702.03542.pdf
"""
bead_magnetic_moment = (((np.pi)*(bead_diameter**3))/(2*perm0))*((perm_bead-perm0)/(perm_bead+(2*perm0)))*Magnetic_Flux_Density  #magnetic moment of the microrobot (volume * Magnitization)

Magnetic_Force = bead_magnetic_moment * Magnetic_Flux_Density                                                                 #This is wrong do not know how to calculation u dot del in python. probably isnt neccesary though

print("Magnetic Force exerted on Microbead =",Magnetic_Force*(1/pico), "pN\n")
# %%
#
#
#
#
#
#
#
#
#
#
#VISCOUS FORCE
# %%
Viscous_Force = (6*np.pi*Viscosity*(bead_diameter/2)*Velocity) * (1/pico)
print("Viscous Force = ",Viscous_Force, "pN")
#%%