#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 14 15:20:47 2021

@author: bizzaro
"""

import numpy as np

milli = (10**(-3))
micro = (10**(-6))
nano = (10**(-9))
pico = (10**(-12))


Pole_Pole_Distance = 980 * micro 
    
L = (Pole_Pole_Distance/2)  
    
measurement_location = [0*micro,0*micro] 







##################################################################################################################
# %%
#VECTOR FORMULATION
 
Left_Pole_Origin = [-L,0]
Top_Pole_Origin = [0,L]
Right_Pole_Origin = [L,0]
Bottom_Pole_Origin = [0,-L]

def r_vector(Pole_Location, measurement_location):
    '''
    returns a 2D array fector in direction of microbot from desired pole
    '''
    x1 = Pole_Location[0]
    y1 = Pole_Location[1]
    
    x2 = measurement_location[0]
    y2 = measurement_location[1]
    r_vector = np.array([(x2-x1),(y2-y1)])
 
    return r_vector

def r_distance(Pole_Location, measurement_location):
    '''
    returns a scalar distance desired pole to microbot
    not the normalized sistance r/L
    '''
    x1 = Pole_Location[0]
    y1 = Pole_Location[1]
    
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
Left_Pole_Vector = r_vector(Left_Pole_Origin, measurement_location) 
Top_Pole_Vector = r_vector(Top_Pole_Origin, measurement_location)
Right_Pole_Vector = r_vector(Right_Pole_Origin, measurement_location)
Bottom_Pole_Vector = r_vector(Bottom_Pole_Origin, measurement_location)
##################################################################################################################
# %%






##################################################################################################################
# %%
#L MATRIX FORMULATION
"""
L matrix calculation. Backbone of force directionality

Bascially its the gradient of the dot product between each pole-bot vector sum(k x4) sum(i x 4) of (ui dot uk ) 


NEED TO ADD FUNCTION THAT DOES THIS FOR EACH VECTOR: r_vect_j_hat = r_vect_j/r_dist

"""
L_matrix = []

r_dist = np.array([Left_Pole_Dist,Top_Pole_Dist,Right_Pole_Dist,Bottom_Pole_Dist])

r_vect_j = np.array([Left_Pole_Vector,Top_Pole_Vector,Right_Pole_Vector,Bottom_Pole_Vector])
r_vect_k = np.array([Left_Pole_Vector,Top_Pole_Vector,Right_Pole_Vector,Bottom_Pole_Vector])

#r_vect_j_hat = r_vect_j/r_dist

for j in range(len(r_vect_j)):
    for k in range(len(r_vect_k)):
        
        
        
        r_dist_hat_j = (np.linalg.norm(r_vect_j[j]))/L
        r_dist_hat_k = (np.linalg.norm(r_vect_k[k]))/L
        
        coef = 1/((r_dist_hat_j**3)*(r_dist_hat_k**3))
        
        
        '''
        replace Left_Pole_Dist with r_dist[i]
        99.9% sure this is correct
        '''
        r_vect_hat_j = (r_vect_j[j]/r_dist[j])
        r_vect_hat_k = (r_vect_k[k]/r_dist[k])
        
        dot_left = (3*(np.dot(r_vect_hat_j  , r_vect_hat_k))
        dot_right = (3*(np.dot(r_vect_hat_j , r_vect_hat_k))
        
     
        Leftside = (1- dot_left/(r_dist_hat_k**2)))*r_vect_hat_k 
        Rightside = (1- dot_right/(r_dist_hat_j**2)))*r_vect_hat_j 
        
        print("Leftside ", Leftside)
        print("Rightside ", Rightside)
        
        Lindex = coef * (Leftside + Rightside)
        
        L_matrix.append(Lindex)
        
        
L_matrix_nump = np.array(L_matrix).reshape((4,4,2))

L_matrix_nump_X = L_matrix_nump[:,:,0]
L_matrix_nump_Y = L_matrix_nump[:,:,1]

print(L_matrix_nump)
##################################################################################################################
# %%


