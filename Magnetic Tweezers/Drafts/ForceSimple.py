#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 13 10:14:32 2021

@author: bizzaro
"""

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"]=(8, 8)

milli = (10**(-3))
micro = (10**(-6))
nano = (10**(-9))
pico = (10**(-12))





"""
force range (∼1–100 pN) from paper.
"""









##################################################################################################################
# %%
#System Properties
perm0 = ((4*np.pi)*(10**(-7)))


Pole_Pole_Distance = 490 * 2 * micro
km = (perm0/4*np.pi)
L = (Pole_Pole_Distance/2) # length to center of workspace

I_max = .9 #A


ILeft = 0 #Left
ITop = 0 #Top
IRight = 1 #Right
IBottom = 0 #Bottom

#Microsphere Properites
Bead_Diameter = (2.8 *micro)
perm_bead = .0000375
 
measurement_location = [0*micro,0*micro]

Volume = ((4/3)*np.pi*((Bead_Diameter/2)**3))
NTurns = 21
reluctance = (1.8*(10**10))
##################################################################################################################
# %%













'''

As the magnetic probe is much more magnetic permeable than air, (μ − μ0 )/(μ + 2μ0 )

'''
##################################################################################################################
# %%

Perm_coef = 1#((perm_bead - perm0)/(perm_bead + (2*perm0)))
kQ = ((3*Volume*(km**2))/(2*perm0*(L**5))) * Perm_coef #coefficient related to the properties of the magnetic particle and themagnetic tweezers
ki = (kQ * ((NTurns/(perm0*reluctance))**2)) #lumped coefficent related to the particle and magnetic curcuit properites
ki_hat = (ki * (I_max**2))
print("Ki = ", ki)
print("Ki_hat = ", ki_hat*(1/pico),"pN\n\n")
##################################################################################################################
#%%


















##################################################################################################################
# %%



Left_Pole_Origin = [-Pole_Pole_Distance/2,0]
Top_Pole_Origin = [0,Pole_Pole_Distance/2]
Right_Pole_Origin = [Pole_Pole_Distance/2,0]
Bottom_Pole_Origin = [0,-Pole_Pole_Distance/2]

def r_vector(pole_pos, measurement_location):
    '''
    returns a 2D array fector in direction of microbot from desired pole
    '''
    x1 = pole_pos[0]
    y1 = pole_pos[1]
    
    x2 = measurement_location[0]
    y2 = measurement_location[1]
    r_vector = np.array([(x2-x1),(y2-y1)])
 
    return r_vector

def r_distance(pole_pos, measurement_location):
    '''
    returns a scalar distance desired pole to microbot
    not the normalized sistance r/L
    '''
    x1 = pole_pos[0]
    y1 = pole_pos[1]
    
    x2 = measurement_location[0]
    y2 = measurement_location[1]
    
    r_dist = (np.sqrt((x2-x1)**2 + (y2-y1)**2))
    return r_dist

#Assigns disances to microbot from each pole to a variable
Left_Pole_Dist = r_distance(Left_Pole_Origin, measurement_location)
Top_Pole_Dist = r_distance(Top_Pole_Origin, measurement_location)
Right_Pole_Dist = r_distance(Right_Pole_Origin, measurement_location)
Bottom_Pole_Dist = r_distance(Bottom_Pole_Origin, measurement_location)

#Assigns direction to microbot from each pole to a variable
Left_Pole_Vector = r_vector(Left_Pole_Origin, measurement_location)/Left_Pole_Dist 
Top_Pole_Vector = r_vector(Top_Pole_Origin, measurement_location)/Top_Pole_Dist
Right_Pole_Vector = r_vector(Right_Pole_Origin, measurement_location)/Right_Pole_Dist
Bottom_Pole_Vector = r_vector(Bottom_Pole_Origin, measurement_location)/Bottom_Pole_Dist
##################################################################################################################
# %%















##################################################################################################################
# %%
#L MATRIX FORMULATION
"""
L matrix calculation. Backbone of force directionality

Bascially its the gradient of the dot product between each pole-bot vector sum(k x4) sum(i x 4) of (ui dot uk ) 
"""
L_matrix = []

#this is r_hat: the nomralized distance
r_dist = [Left_Pole_Dist/L,Top_Pole_Dist/L,Right_Pole_Dist/L,Bottom_Pole_Dist/L]

r_vect_j = [Left_Pole_Vector,Top_Pole_Vector,Right_Pole_Vector,Bottom_Pole_Vector]
r_vect_k = [Left_Pole_Vector,Top_Pole_Vector,Right_Pole_Vector,Bottom_Pole_Vector]


for j in r_vect_j:
    for k in r_vect_k:
        
        rdj = (np.linalg.norm(j))
        rdk = (np.linalg.norm(k))
        
        coef = 1/((rdj**3)*(rdk**3))
        
        Leftside = (1- (3*(np.dot(j,k))/(rdk**2)))*k
        Rightside = (1- (3*(np.dot(j,k))/(rdj**2)))*j
        
        Lindex = coef * (Leftside + Rightside)
        
        L_matrix.append(Lindex)
        
        
L_matrix_nump = np.array(L_matrix).reshape((4,4,2))

L_matrix_nump_X = L_matrix_nump[:,:,0]
L_matrix_nump_Y = L_matrix_nump[:,:,1]



##################################################################################################################
# %%















##################################################################################################################
# %%
#N MATRIX FORMULATION 
'''
N_matrix calcualtions. N = KI_T  * L  *  KI
KI matrix --> charge distributuons throughtout magnetic curcuit
'''
KI = np.array([[3/4, -1/4, -1/4, -1/4],  #matrix KI is the distribution matrix of the magnetic flux.
              [-1/4, 3/4, -1/4, -1/4],
              [-1/4, -1/4, 3/4, -1/4],
              [-1/4, -1/4, -1/4, 3/4]])


N_Matrix_X = np.matmul(np.matmul(KI.transpose() , L_matrix_nump_X), KI)
N_Matrix_Y = np.matmul(np.matmul(KI.transpose() , L_matrix_nump_Y), KI)
##################################################################################################################
#%%












'''
F-matrix Calcualtions. F = ki * Imat * Nmat * Imat   =   ki * Imax^2  * Ihatmat * N * Ihat = kihat * Fhat
'''
##################################################################################################################
# %%
#this force represents teh complete force model but does not consider magnetic pole saturation 
print("This should matches the force on microbot from 'MagneticForceAnalysis.py' code which directly uses Q_matrix")
I_matrix = np.array([ILeft,ITop,IRight,IBottom])

Fx = ki * np.matmul(np.matmul(I_matrix, N_Matrix_X), I_matrix.transpose())
Fy = ki * np.matmul(np.matmul(I_matrix, N_Matrix_Y), I_matrix.transpose())

print("Fx = ", Fx*(1/pico), "pN")
print("Fy = ", Fy*(1/pico), "pN\n\n")
##################################################################################################################
# %%









#still have issue of normalized force matching regular force







##################################################################################################################
#%%
# F_HAT MATRIX FORMULATION
'''
F_hat Calcuations. Normalized force vector which becomes dimensionless.
    As magnetic poles saturate when excessive current is applied,
    the input current is therefore limited to I_max  , which is 
    the maximum current. Normalizing the current vector by dividing 
    I_matrix by I_max , the force model can be rewritten as the following
'''

I_matrix_hat = I_matrix/I_max

F_hat_X = np.matmul(np.matmul(I_matrix_hat, N_Matrix_X), I_matrix_hat.transpose())
F_hat_Y = np.matmul(np.matmul(I_matrix_hat, N_Matrix_Y), I_matrix_hat.transpose())

print("Actual force nomalized by maximum currnet")

FORCEX = ki_hat * F_hat_X
FORCEY = ki_hat * F_hat_Y

print("Fx = ", FORCEX*(1/pico), "pN")
print("Fy = ", FORCEY*(1/pico), "pN\n")

print("Fx_hat_direction = ",F_hat_X)
print("Fy_hat_direction = ",F_hat_Y, "\n")


##################################################################################################################
#%%































'''
Plot visual model
'''
##################################################################################################################
# %%
#plot microbot
plt.scatter(measurement_location[0],measurement_location[1],s=100,facecolors = "k")
 
#plot the vectors from each pole to microbot    
plt.quiver(Left_Pole_Origin[0],Left_Pole_Origin[1] ,Left_Pole_Vector[0] , Left_Pole_Vector[1],color='r', scale = None)
plt.quiver(Top_Pole_Origin[0],Top_Pole_Origin[1] ,Top_Pole_Vector[0] , Top_Pole_Vector[1],color='b',scale = None)
plt.quiver(Right_Pole_Origin[0],Right_Pole_Origin[1] ,Right_Pole_Vector[0] , Right_Pole_Vector[1],color='g',scale = None)
plt.quiver(Bottom_Pole_Origin[0],Bottom_Pole_Origin[1] ,Bottom_Pole_Vector[0] , Bottom_Pole_Vector[1],color='y',scale = None)

#Plot force direction on microbot
plt.quiver(measurement_location[0],measurement_location[1] ,Fx, Fy,color='k', scale = None)
plt.quiver(measurement_location[0],measurement_location[1] ,FORCEX, FORCEY,color='b', scale = 2)


#plt.quiver(Left_Pole_Origin, V[0], V[1], color=['r','b','g'])
plt.xlim(-(Pole_Pole_Distance/2),(Pole_Pole_Distance/2))
plt.ylim(-(Pole_Pole_Distance/2),(Pole_Pole_Distance/2))
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.show()
##################################################################################################################
# %%













'''
Quick force simplification check from paper

     Accord- ing to the definition of magnetic charges, the summation of the 
     magnetic charges is thus equal to zero, q1 +q2 +q3 +q4 =0.  In other words, 
     no excessive magnetic charge exists in the system. This singularity of KI 
     indicates that there are multiple solutions associated with the same Q. 
     If a similar constraint, I1 +I2 +I3 +I4 =0 is applied to the input currents, 
     a one-to-one relationship between Q and I can be established. The force model 
     can then be simplified. In particular, at the center of the workspace the 
     force model is greatly abridged due to the symmetry of the setup.
'''
##################################################################################################################
#%%
print("##################################################################################################################\n")
print("If  I1 + I2 + I3 + I4 = 0, then the force on the microbot at model at the workspace center is the following\n")
#only if I1+I2+I3+I4 =0
Fx_hat_0_0 = 6*((I_matrix_hat[0]**2) - (I_matrix_hat[2]**2))
Fy_hat_0_0 = 6*((I_matrix_hat[1]**2) - (I_matrix_hat[3]**2))


print("Fx_hat_direction at (0,0) = ", Fx_hat_0_0)
print("Fy_hat_direction at (0,0) = ", Fy_hat_0_0, "\n")

#Therefore, the range of the magnetic force at the center of the workspace is [ −6ki_hatˆ 6ki_hat ] in each direction since |I_hat | ≤ 1.
print("And the range of magnetic force at the center is\n")
Fx_mag_0_0 = -6 * ki_hat
Fy_mag_0_0 =  6 * ki_hat

print("Fx_magnitude at (0,0) = ", Fx_mag_0_0*(1/pico),"pN")
print("Fy_magnitude at (0,0) = ", Fy_mag_0_0*(1/pico),"pN")
##################################################################################################################
#%%