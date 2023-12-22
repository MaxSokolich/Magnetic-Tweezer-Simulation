import numpy as np
import matplotlib.pyplot as plt




#VECTOR FORMULATION
##################################################################################################################
# %%
milli = (10**-3)
Pole_Pole_Distance = 2 * milli

microrobot_location = (.5*milli,.5*milli)
Left_Pole_Origin = [-Pole_Pole_Distance/2,0]
Top_Pole_Origin = [0,Pole_Pole_Distance/2]
Right_Pole_Origin = [Pole_Pole_Distance/2,0]
Bottom_Pole_Origin = [0,-Pole_Pole_Distance/2]

def r_vector(Pole_Location, microrobot_location):
    '''
    returns a 2D array fector in direction of microbot from desired pole
    '''
    x1 = Pole_Location[0]
    y1 = Pole_Location[1]
    
    x2 = microrobot_location[0]
    y2 = microrobot_location[1]
    r_vector = np.array([(x2-x1),(y2-y1)])
 
    return r_vector

def r_distance(Pole_Location, microrobot_location):
    '''
    returns a scalar distance desired pole to microbot
    not the normalized sistance r/L
    '''
    x1 = Pole_Location[0]
    y1 = Pole_Location[1]
    
    x2 = microrobot_location[0]
    y2 = microrobot_location[1]
    
    r_dist = (np.sqrt((x2-x1)**2 + (y2-y1)**2))
    return r_dist

#Assigns disances to microbot from each pole to a variable
Left_Pole_Dist = r_distance(Left_Pole_Origin, microrobot_location)
Top_Pole_Dist = r_distance(Top_Pole_Origin, microrobot_location)
Right_Pole_Dist = r_distance(Right_Pole_Origin, microrobot_location)
Bottom_Pole_Dist = r_distance(Bottom_Pole_Origin, microrobot_location)

#Assigns direction to microbot from each pole to a variable
Left_Pole_Vector = r_vector(Left_Pole_Origin, microrobot_location) 
Top_Pole_Vector = r_vector(Top_Pole_Origin, microrobot_location)
Right_Pole_Vector = r_vector(Right_Pole_Origin, microrobot_location)
Bottom_Pole_Vector = r_vector(Bottom_Pole_Origin, microrobot_location)

#plot microbot
plt.figure()
plt.scatter(microrobot_location[0],microrobot_location[1],s=100,facecolors = "k")
 
#plot the vectors from each pole to microbot    
plt.quiver(Left_Pole_Origin[0],Left_Pole_Origin[1] ,Left_Pole_Vector[0] , Left_Pole_Vector[1],color='r', scale = None)
plt.quiver(Top_Pole_Origin[0],Top_Pole_Origin[1] ,Top_Pole_Vector[0] , Top_Pole_Vector[1],color='b',scale = None)
plt.quiver(Right_Pole_Origin[0],Right_Pole_Origin[1] ,Right_Pole_Vector[0] , Right_Pole_Vector[1],color='g',scale = None)
plt.quiver(Bottom_Pole_Origin[0],Bottom_Pole_Origin[1] ,Bottom_Pole_Vector[0] , Bottom_Pole_Vector[1],color='y',scale = None)
##################################################################################################################
# %%




#Magnetic field matrix from 4 poles
const = (4 *np.pi) / (4*np.pi*10**(-7))
#B1 =  const * (Q1/(r1**2)) 
