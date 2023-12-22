milli = (10**(-3))
micro = (10**(-6))
nano = (10**(-9))
pico = (10**(-12))


Pole_Pole_Distance = 1000 * micro                # pole to pole distance in um, defines workspace area
measurement_location = [100*micro,0*micro]         # microbot position of interest
workspace_area = [-900*micro, 900*micro]


IBottom = 0                                          # Applied current to left coil/pole in Amps
ITop = 0                                         # Applied current to top coil/pole in Amps
IRight = 1                                      # Applied.5urrent to right coil/pole in Amps
ILeft= 0                                      # Applied current to bottom coil/pole in Amps          
        

#OUTPUTS

#FORCE ON MICROBOT, MAGNETIC FIELD , VELOCITY OF MICROBOT AT MEASURMENT LOCATION and PRESCIRBED INPUT CURRENTS
#GRAPH OF MANGETIC FIELD LINES AND MICROBOT AT MEASURMENT LOCATION
#GRAPH OF FORCE FIELD FOR PRESCRIBED CURRENT COMBINATION
#GRAPH OF FORCE ENVELOPE ON MICROBOT AT PRESCRIBED MEASUREMENT LOCATION

#TO DO
'''
-L direction matrix still buggy. outputing force directions that dont make sense for complex current combinations
-Force envelope vector magnetiudes. cant figure out how to scale the lengths

-Convert to 3D
-Convert to Current inputs
'''




#DESCRIPTION
'''


MAGNETIC FORCE ANALYSIS / MAGNETIC MONOPOLE APPROXIMATION / LUMPED-ELEMENT MODEL

MAGNETIC MONOPOLE APPRIXMATION THEORY FROM: 
    "Design, Implementation, and Force Modeling of Quadrupole Magnetic Tweezers", 
    
    [1] -  https://ieeexplore.ieee.org/document/5280242
    [2] -  https://sites.google.com/site/zpzhang/research-highlights/magnetic-tweezers

    The lumped-element model (also called lumped-parameter model, or 
    lumped-component model) simplifies the description of the 
    behaviour of spatially distributed physical systems into
    a topology consisting of discrete entities that approximate 
    the behaviour of the distributed system under certain assumptions.
    It is useful in electrical systems (including electronics), 
    mechanical multibody systems, heat transfer, acoustics, etc.

    When viewing magnetic poles of the tweezers from the workspace, 
    they are sharp tips and generate a magnetic field that looks to 
    the magnetic particle as though it is generated by point sources. 
    Magnetic monopole approximation [25], [30] is there- fore employed 
    in this analysis to model the magnetic field gen- erated by the 
    quadrupole magnetic tweezers in the workspace. In this model, 
    the magnetic field generated by each magnetic pole is approximated 
    by the field of a point magnetic charge associated with the magnetic 
    pole, and the total magnetic field produced by the system is obtained 
    by applying the principle of superposition. Although a magnetic monopole 
    is a hypothetical point source, this approximation greatly simplifies the 
    process of modeling and yields an expression of the magnetic field in 
    the workspace similar to that of an electric field.More importantly,this 
    expression enables us later to derive a compact analytical force
    model that accurately characterizes the nonlinearity of the magnetic 
    force exerting on the magnetic particle with respect to the applied 
    currents to the coils and the position dependency of the magnetic 
    force in the workspace.

    When a magnetic particle is placed inside a magnetic field, 
    it will be magnetized, and in turn alter its surrounding magnetic field. 
    The magnetic force experienced by the magnetic particle is then determined 
    by the particle’s magnetization and the altered field. To calculate the 
    magnetic force using the original magnetic field, the effective magnetization 
    of the particle should be used. For a spherical par- ticle, 
    its effective magnetization is proportional to the external magnetic field [31] 
    until the magnetization saturates.
'''
#INITLIZIATION

#IMPORTS
# %%
import numpy as np
import matplotlib.pyplot as plt
import random as random
from tqdm import tqdm
import matplotlib.colors as colors

plt.rcParams["figure.figsize"]=(8, 8)


perm0 = ((4*np.pi)*(10**(-7)))
# %%




#INPUTS TO SYSTEM
##################################################################################################################
#%%


                                  

#System Properties
                 
L = (Pole_Pole_Distance/2)                         # length to center of workspace
Gradient_Space = (L/4)                             # Maximum gradient estimation distance 
I_max = 1.5                                         # Maximum current to avoid pole saturation --> value taken from paper
NTurns = 400                                        # number of turns on a single coil --> taken from paper
#reluctance = (1.8*(10**9))                         # magnetic resitance of air gap --> value from paper 



#MICROSPHERE PROPERTIES
Bead_Diameter = (23.7 * micro)                       # bead diameter in micrometers
perm_bead = .000375                               # neglected for force calculations: see "magnetic coeffceitns for why"
Volume = ((4/3)*np.pi*((Bead_Diameter/2)**3))      # volume of microsphere



#This wont really change anything
Air_AreaC = (.18*milli*milli)                     # cross sectioal area of air gap 



gridp = 50                                         # grid resolution of magnetic streamlines

#num_arrows = 1000                                    # Number of force vectors shown on force envelope

#DRAG FORCE INPUTS
viscosity = (8.90 * (10**(-4)))                    # Water Pa·s
velocity_X = (27)*micro                            # x veclity component of microsphere
velocity_Y = (0)*micro                             # y velocity component of micrpsphere
##################################################################################################################
# %%










#MAGNETIC CURCUIT ANALYSIS
##################################################################################################################
# %%
#print("##################################################################################################################")
#print("                                        MAGNETIC CURCUIT ANALYSIS\n")

'''
In this circuit, Ra is significantly larger than Rp and Ry 
as the yoke and poles are made of materials much more magnetic permeable than air.
'''
#Yoke Reluctance
yoke_length = (0 *milli)
YpermR = 5000
Yoke_AreaC = (25*milli*milli)
Ry = (yoke_length)/(YpermR*perm0*Yoke_AreaC)
#print("Yoke Reluctance = ", Ry)



#Pole Reluctance
pole_length = (20 *milli)
PpermR = 740
Pole_AreaC = (1*milli*milli)
Rp = (pole_length)/(PpermR*perm0*Pole_AreaC)
#print("Pole Reluctance = ", Rp)

#Air Gap Reluctance (Length measured from pole to workspace center)
Air_length = (Pole_Pole_Distance/2)
ApermR = 1

'''
#-------IMPORTANT--------the air cross sectionarea is what should be calculted for as              
the relutance of the air gap was calcualted from ANSYS simulation. It is hard to calculate this because
the magnetic field is spread out 

Thus the model will fail if pole diameter is 2 big for 2 reasons;  monopole approximation doesnt hold, 
it will be harder to model the magnetic curcuit and obtain the relationship between input current 
and magnetic flux
'''
Ra = (Air_length)/(ApermR*perm0*Air_AreaC)

reluctance = Ra + Rp


#print("Air Gap Reluctance = ", Ra, "\n\n")
##################################################################################################################
#%%











#MAGNETIC COEFFICIENTS
##################################################################################################################
# %%
#print("##################################################################################################################")
#print(                                        "MAGNETIC COEFFICENTS \n")
'''
ji_hat is basically the force gain and is mulitply my the direction because the magnetic force is dependent on L matrix or where the particle is
As the magnetic probe is much more magnetic permeable than air, (μ − μ0 )/(μ + 2μ0 ) ~ 1

'''

Perm_coef = 1#((perm_bead - perm0)/(perm_bead + (2*perm0)))

km = (perm0/(4*np.pi))
kQ = ((3*Volume*(km**2))/(2*perm0*(L**5))) * Perm_coef #coefficient related to the properties of the magnetic particle and themagnetic tweezers
ki = (kQ * ((NTurns/(perm0*reluctance))**2)) #lumped coefficent related to the particle and magnetic curcuit properites
ki_hat = (ki * (I_max**2))


#print("Ki = ", ki)
print("Ki_hat = ", ki_hat*(1/pico),"pN\n\n")
##################################################################################################################
#%%










#VECTOR FORMULATION
##################################################################################################################
# %%

 
Left_Pole_Origin = [-Pole_Pole_Distance/2,0]
Top_Pole_Origin = [0,Pole_Pole_Distance/2]
Right_Pole_Origin = [Pole_Pole_Distance/2,0]
Bottom_Pole_Origin = [0,-Pole_Pole_Distance/2]

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










#Q MATRIX FORMULATION, 
##################################################################################################################
#%%

"""
From qi = magnetic flux/perm0 

a Q matrix is born which equals = (NTurns/(perm0*reluctance))    *KIMatrix @ Imatrix

NOTE: from Zhang thesis : "As the six poles in our system are not connected to the same yoke as in [71] , the KI matrix in the vector of magnetic charges Q can be neglected."
"""

'''

KI_matrix = np.array([[3/4, -1/4, -1/4, -1/4],  #matrix KI is the distribution matrix of the magnetic flux.
                      [-1/4, 3/4, -1/4, -1/4],
                      [-1/4, -1/4, 3/4, -1/4],
                      [-1/4, -1/4, -1/4, 3/4]])

    
'''
KI_matrix = np.array([[3/4, -1/4, -1/4, -1/4],  #matrix KI is the distribution matrix of the magnetic flux.
                      [-1/4, 3/4, -1/4, -1/4],
                      [-1/4, -1/4, 3/4, -1/4],
                      [-1/4, -1/4, -1/4, 3/4]])

I_matrix = np.array([ILeft,ITop,IRight,IBottom])


Q_Matric_Const = (NTurns/(perm0*reluctance))
Qmatrix = Q_Matric_Const * np.matmul(KI_matrix , I_matrix)
print(Qmatrix)
Q_Left = Qmatrix[0]
Q_Top = Qmatrix[1]
Q_Right = Qmatrix[2]
Q_Bottom = Qmatrix[3]
##################################################################################################################
#%%










#MAGNETIC FIELD ANALYSIS
##################################################################################################################
# %%

print("##################################################################################################################")
print("                                        MAGNETIC FIELD ANALYSIS\n")



def Bfield_X(q,Pole_Location, Meas_Loc):
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


#thru superspoistion principle sum up magnetic field contribution from each pole and measure the field at the location
BTX = (1000)*(Bfield_X(Q_Left,Left_Pole_Origin, measurement_location) + Bfield_X(Q_Top,Top_Pole_Origin, measurement_location) + Bfield_X(Q_Right,Right_Pole_Origin, measurement_location) +  Bfield_X(Q_Bottom,Bottom_Pole_Origin, measurement_location))
BTY = (1000)*(Bfield_Y(Q_Left,Left_Pole_Origin, measurement_location) + Bfield_Y(Q_Top,Top_Pole_Origin, measurement_location) + Bfield_Y(Q_Right,Right_Pole_Origin, measurement_location) +  Bfield_Y(Q_Bottom,Bottom_Pole_Origin, measurement_location))
BT = [BTX,BTY]
print("\nLocation = ",[measurement_location[0]*(1/micro),measurement_location[1]*(1/micro)] , "um")
print("Magnetic Field =", BT, 'mT\n\n')



#STREAMPLOT CONTROUR PLOT
X = np.linspace(-(Pole_Pole_Distance/2), (Pole_Pole_Distance/2), gridp)
Y = np.linspace(-(Pole_Pole_Distance/2),(Pole_Pole_Distance/2), gridp)
x,y = np.meshgrid(X, Y)


#Thru Super Positon

BTotal_X = (1000)*(Bfield_X(Q_Left,Left_Pole_Origin,[x,y]) + Bfield_X(Q_Top,Top_Pole_Origin,[x,y]) + Bfield_X(Q_Right,Right_Pole_Origin,[x,y]) +  Bfield_X(Q_Bottom,Bottom_Pole_Origin,[x,y]))
BTotal_Y = (1000)*(Bfield_Y(Q_Left,Left_Pole_Origin,[x,y]) + Bfield_Y(Q_Top,Top_Pole_Origin,[x,y]) + Bfield_Y(Q_Right,Right_Pole_Origin,[x,y]) +  Bfield_Y(Q_Bottom,Bottom_Pole_Origin,[x,y]))



##################################################################################################################
# %%










#MAGNETIC FIELD VS DISTANCE TO CENTER
##################################################################################################################
#%%



B_Dist_X = []

dist_list = np.linspace(0, L)

for i in dist_list:
    B_dist = (1000)*(Bfield_X(Q_Right, Right_Pole_Origin, [i,0]))
    B_Dist_X.append(B_dist)

beta = (1/micro)
dist_list_inv = np.linspace(L, .000000001)
Inverse_Square = (1/(beta*(dist_list**2)))
##################################################################################################################
#%%











#L MATRIX FORMULATION
##################################################################################################################
# %%

"""
L matrix calculation. Backbone of force directionality

Bascially its the gradient of the dot product between each pole-bot vector sum(k x4) sum(i x 4) of (ui dot uk ) 


NEED TO ADD FUNCTION THAT DOES THIS FOR EACH VECTOR: r_vect_j_hat = r_vect_j/r_dist

"""
L_matrix = []

r_dist = [Left_Pole_Dist,Top_Pole_Dist,Right_Pole_Dist,Bottom_Pole_Dist]

r_vect_j = [Left_Pole_Vector,Top_Pole_Vector,Right_Pole_Vector,Bottom_Pole_Vector]
r_vect_k = [Left_Pole_Vector,Top_Pole_Vector,Right_Pole_Vector,Bottom_Pole_Vector]

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
        
        Leftside = (1- (3*(np.dot(r_vect_hat_j  ,r_vect_hat_k))/(r_dist_hat_k**2)))*r_vect_hat_k 
        Rightside = (1- (3*(np.dot(r_vect_hat_j  ,r_vect_hat_k))/(r_dist_hat_j**2)))*r_vect_hat_j 
        
       
        
        Lindex = coef * (Leftside + Rightside)
        
        L_matrix.append(Lindex)
        
        
L_matrix_nump = np.array(L_matrix).reshape((4,4,2))

L_matrix_nump_X = L_matrix_nump[:,:,0]
L_matrix_nump_Y = L_matrix_nump[:,:,1]

##################################################################################################################
# %%










#N MATRIX FORMULATION 
##################################################################################################################
# %%

'''
N_matrix calcualtions. N = KI_T  * L  *  KI
KI matrix --> charge distributuons throughtout magnetic curcuit
'''

N_Matrix_X = np.matmul(np.matmul(KI_matrix.transpose() , L_matrix_nump_X), KI_matrix)
N_Matrix_Y = np.matmul(np.matmul(KI_matrix.transpose() , L_matrix_nump_Y), KI_matrix)

#print([N_Matrix_X,N_Matrix_Y])
##################################################################################################################
#%%










#GLOBAL F_MAGNETIC VECTOR FORMULATION
##################################################################################################################
# %%

print("##################################################################################################################")
print("                                        F PLAIN MATRIX ANALYSIS\n")
'''
F-matrix Calcualtions. F = ki * Imat * Nmat * Imat   =   ki * Imax^2  * Ihatmat * N * Ihat = kihat * Fhat
'''

FQx = np.matmul(np.matmul(Qmatrix,   L_matrix_nump_X) , Qmatrix.transpose())
FQy = np.matmul(np.matmul(Qmatrix,   L_matrix_nump_Y) , Qmatrix.transpose())

Fqx  = FQx * kQ
Fqy =  FQy * kQ





#this force represents teh complete force model but does not consider magnetic pole saturation 
I_matrix = np.array([ILeft,ITop,IRight,IBottom])

Fx = ki * np.matmul(np.matmul(I_matrix, N_Matrix_X), I_matrix.transpose())
Fy = ki * np.matmul(np.matmul(I_matrix, N_Matrix_Y), I_matrix.transpose())

print("F_MAG_X = ", Fx*(1/pico), "pN")
print("F_MAG_Y = ", Fy*(1/pico), "pN\n\n")
##################################################################################################################
# %%










# DRAG FORCE CALCULATION
##################################################################################################################
# %%

#Force due to stokes drag on a sphere
F_DRAG_X = 3*np.pi*viscosity*(Bead_Diameter)*velocity_X * (1/pico)
F_DRAG_Y = 3*np.pi*viscosity*(Bead_Diameter)*velocity_Y * (1/pico)
print("DRAG FORCE ON MICROSPHERE FOR VELOCITY = ", [velocity_X*(1/micro), velocity_Y*(1/micro)], "um/s")
print("F_DRAG_X = ",F_DRAG_X, "pN")
print("F_DRAG_Y = ",F_DRAG_Y, "pN\n")
##################################################################################################################
# %%










#VELOCITY CALCULATION
##################################################################################################################
#%%

'''
F_magnetic => F_drag to induce motion 
If F_drag = 3*np.pi*viscosity*(Bead_Diameter)*velocity_X * (1/pico)

then V = F_magnetic / (3*np.pi*viscosity*(Bead_Diameter))
'''
VEL_X = (Fx/(3*np.pi*viscosity*(Bead_Diameter)))
VEL_Y = (Fy/(3*np.pi*viscosity*(Bead_Diameter)))


print("VELOCITY CALCULATION FROM GIVEN FORCE")
print("VELOCITY X MAG = ",VEL_X*(1/micro), "um/s")
print("VELOCITY Y MAG = ",VEL_Y*(1/micro), "um/s\n")

##################################################################################################################
#%%











# F_HAT VECTOR FORMULATION
##################################################################################################################
#%%

#print("##################################################################################################################")
#print("                                        F HAT MATRIX ANALYSIS\n")
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

#print("Force nomalized by maximum currnet")



#print("Fx_hat = ", F_hat_X*(1/pico), "pN")
#print("Fy_hat = ", F_hat_Y*(1/pico), "pN\n")

print("Fx_hat_direction = ",ki_hat*F_hat_X)
#print("Fy_hat_direction = ",F_hat_Y, "\n\n")


##################################################################################################################
#%%










#PLOTTING
##################################################################################################################
# %%


#plot microbot
plt.figure()
plt.scatter(measurement_location[0],measurement_location[1],s=100,facecolors = "k")
 
#plot the vectors from each pole to microbot    
plt.quiver(Left_Pole_Origin[0],Left_Pole_Origin[1] ,Left_Pole_Vector[0] , Left_Pole_Vector[1],color='r', scale = None)
plt.quiver(Top_Pole_Origin[0],Top_Pole_Origin[1] ,Top_Pole_Vector[0] , Top_Pole_Vector[1],color='b',scale = None)
plt.quiver(Right_Pole_Origin[0],Right_Pole_Origin[1] ,Right_Pole_Vector[0] , Right_Pole_Vector[1],color='g',scale = None)
plt.quiver(Bottom_Pole_Origin[0],Bottom_Pole_Origin[1] ,Bottom_Pole_Vector[0] , Bottom_Pole_Vector[1],color='y',scale = None)

#Plot force direction on microbot
plt.quiver(measurement_location[0],measurement_location[1] ,F_hat_X, F_hat_Y,color='k', scale = None)

#plot magnetic field streamlines
speed = np.sqrt((BTotal_X)**2+(BTotal_Y)**2)
plt.streamplot(x,y,BTotal_X,BTotal_Y, color=speed,cmap='coolwarm',density = 2,norm = colors.LogNorm(vmin=speed.min(), vmax=speed.max() ))
#plt.colorbar()


#uncomment for contour plot
'''
cs = plt.contourf(x,y,speed,cmap='coolwarm')
plt.colorbar()
'''

plt.xlim(-(Pole_Pole_Distance/2),(Pole_Pole_Distance/2))
plt.ylim(-(Pole_Pole_Distance/2),(Pole_Pole_Distance/2))
plt.title("Magnetic Field Lines and Magnetic Force Direction on Microbot")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.show()






#plot magnetic field vs distance from pole to workspace center
plt.figure()
plt.plot(dist_list,B_Dist_X, color = "r", linewidth = 5,label = "Magnetic Field (mT)")
plt.plot(dist_list_inv ,Inverse_Square, color = "b", linewidth = 5, label = ("Inverse Square Law, beta = %s"%(beta)))

plt.legend()
plt.xlim(0,(Pole_Pole_Distance/2))
plt.ylim(0,(max(B_Dist_X)))

plt.title("Magnetic Field Strength (mt) vs Distance from Right Pole")
plt.xlabel("X (m) (where 0 indicates workspace center)")
plt.ylabel("B (mT)")

plt.show()

##################################################################################################################
#%%










#FORCE SIMPLIFICATION
"""
##################################################################################################################
#%%
#FORCE SIMPLIFICATION
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
"""










#FORCE FIELD
##################################################################################################################
#%%
F_Field_X = []
F_Field_Y = []

rob_pos = []
for i in (range(len(x))):
    for j in range(len(y)):
      
        
        #Assigns disances to microbot from each pole to a variable
        
        robot_x_pos = x[i][j]
        robot_y_pos = y[i][j]
        
        #print(robot_x_pos)
        robot_pos = [robot_x_pos,robot_y_pos]
        
        Left_Pole_Dist = r_distance(Left_Pole_Origin, robot_pos)
        Top_Pole_Dist = r_distance(Top_Pole_Origin, robot_pos)
        Right_Pole_Dist = r_distance(Right_Pole_Origin, robot_pos)
        Bottom_Pole_Dist = r_distance(Bottom_Pole_Origin, robot_pos)
        
        #Assigns direction to microbot from each pole to a variable
        Left_Pole_Vector = r_vector(Left_Pole_Origin, robot_pos) 
        Top_Pole_Vector = r_vector(Top_Pole_Origin, robot_pos)
        Right_Pole_Vector = r_vector(Right_Pole_Origin, robot_pos)
        Bottom_Pole_Vector = r_vector(Bottom_Pole_Origin, robot_pos)
         
        #L MATRIX FORMULATION
        '''
        L matrix calculation. Backbone of force directionality
        
        Bascially its the gradient of the dot product between each pole-bot vector sum(k x4) sum(i x 4) of (ui dot uk ) 
        
        
        NEED TO ADD FUNCTION THAT DOES THIS FOR EACH VECTOR: r_vect_j_hat = r_vect_j/r_dist
        
        '''
        L_matrix = []
        
        r_dist = [Left_Pole_Dist,Top_Pole_Dist,Right_Pole_Dist,Bottom_Pole_Dist]
        
        r_vect_j = [Left_Pole_Vector,Top_Pole_Vector,Right_Pole_Vector,Bottom_Pole_Vector]
        r_vect_k = [Left_Pole_Vector,Top_Pole_Vector,Right_Pole_Vector,Bottom_Pole_Vector]
        
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
                
                Leftside = (1- (3*(np.dot(r_vect_hat_j  ,r_vect_hat_k))/(r_dist_hat_k**2)))*r_vect_hat_k 
                Rightside = (1- (3*(np.dot(r_vect_hat_j  ,r_vect_hat_k))/(r_dist_hat_j**2)))*r_vect_hat_j 
                
               
                
                Lindex = coef * (Leftside + Rightside)
                
                L_matrix.append(Lindex)
                
                
        L_matrix_nump = np.array(L_matrix).reshape((4,4,2))
        
        L_matrix_nump_X = L_matrix_nump[:,:,0]
        L_matrix_nump_Y = L_matrix_nump[:,:,1]
          
        
        
      
        
        N_Matrix_X = np.matmul(np.matmul(KI_matrix.transpose() , L_matrix_nump_X), KI_matrix)
        N_Matrix_Y = np.matmul(np.matmul(KI_matrix.transpose() , L_matrix_nump_Y), KI_matrix)
        
           
        I_matrix = np.array([ILeft,ITop,IRight,IBottom])
        
        Fx = ki * np.matmul(np.matmul(I_matrix, N_Matrix_X), I_matrix.transpose())
        Fy = ki * np.matmul(np.matmul(I_matrix, N_Matrix_Y), I_matrix.transpose())
        F_Field_X.append(Fx*(1/pico))
        F_Field_Y.append(Fy*(1/pico))
        
   
        
F_Field_X_arr = np.array(F_Field_X).reshape(gridp,gridp)
F_Field_Y_arr = np.array(F_Field_Y).reshape(gridp,gridp)



plt.figure()
Fspeed = np.sqrt((F_Field_X_arr)**2+(F_Field_Y_arr)**2)
plt.streamplot(x,y,F_Field_X_arr, F_Field_Y_arr, color=Fspeed,cmap='coolwarm',density = 2, norm = colors.LogNorm(vmin=Fspeed.min(), vmax=Fspeed.max() ))


   
#force field x direction
  
plt.title("FORCE FIELD OF WORKSPACE: 1 Amp Top Pole")
plt.xlim(-Pole_Pole_Distance/2,(Pole_Pole_Distance)/2)
plt.ylim(-Pole_Pole_Distance/2,(Pole_Pole_Distance)/2)
plt.xlabel("x (um)")
plt.ylabel("y (um)")
plt.scatter(measurement_location[0],measurement_location[1],s=100,facecolors = "k")
plt.quiver(measurement_location[0],measurement_location[1] ,F_hat_X, F_hat_Y,color='k', scale = None)
plt.show()

'''
F_Dist = []
x = np.linspace(workspace_area[0],workspace_area[1],gridp)

for i in range(len(x)):
    term1 = F_Field_X_arr[24][i]
    term2 = F_Field_Y_arr[24][i]
    
    FMAG = np.sqrt(term1**2 + term2**2)
  

    F_Dist.append(FMAG)


plt.figure()
plt.plot(x,F_Dist, color = "r", linewidth = 5)


plt.xlim(workspace_area[0],workspace_area[1])
plt.ylim(0,100)



plt.title("Force Field vs Distance from Right Pole")
plt.xlabel("X (m) (where 0 indicates workspace center)")
plt.ylabel("F (pN)")
plt.yscale("log")
plt.show()
'''


##################################################################################################################
#%%












#FORCE ENVELOPE
##################################################################################################################
#%%

'''
To characterize the anisotropy of force generation, 
a measure Γ is defined as the ratio between the smallest 
magnitude and the largest magnitude on the force envelope. 
It is evident that Γ has the largest value, 0.7, at the center 
point and decreases while approaching each of the four poles. 
Fig. 7 shows the contour of Γ. It can be seen that Γ decreases 
rapidly when the magnetic particle moves from the center toward a
 magnetic pole. When Γ is too small, it may not be practically 
 possible to generate the desired force in the weakest force 
 direction due to various uncertainties of the system, such as 
 assembly errors of the tweezers and fluctuations of the commanded 
 currents. It is therefore important to limit the workspace to be within the 
 region where Γ > 0.1 for effective manipulation.
'''
FXX = []
FYY = []


num_arrows = 1000
current2 = np.linspace(0,1,num_arrows)
arrow_start_X = np.full(num_arrows,measurement_location[0]*(1/micro))
arrow_start_Y = np.full(num_arrows,measurement_location[1]*(1/micro))


for j in (range(len(current2))):
      

    ILeft = random.choice(current2)                                           # Applied current to left coil/pole in Amps
    ITop = random.choice(current2)                                         # Applied current to top coil/pole in Amps
    IRight = random.choice(current2)                                        # Applied current to right coil/pole in Amps
    IBottom = random.choice(current2)                                         # Applied current to bottom coil/pole in Amps          
                   
 
    
    
    
    
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
 
    #L MATRIX FORMULATION
    '''
    L matrix calculation. Backbone of force directionality
    
    Bascially its the gradient of the dot product between each pole-bot vector sum(k x4) sum(i x 4) of (ui dot uk ) 
    
    
    NEED TO ADD FUNCTION THAT DOES THIS FOR EACH VECTOR: r_vect_j_hat = r_vect_j/r_dist
    
    '''
    L_matrix = []
    
    r_dist = [Left_Pole_Dist,Top_Pole_Dist,Right_Pole_Dist,Bottom_Pole_Dist]
    
    r_vect_j = [Left_Pole_Vector,Top_Pole_Vector,Right_Pole_Vector,Bottom_Pole_Vector]
    r_vect_k = [Left_Pole_Vector,Top_Pole_Vector,Right_Pole_Vector,Bottom_Pole_Vector]
    
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
            
            Leftside = (1- (3*(np.dot(r_vect_hat_j  ,r_vect_hat_k))/(r_dist_hat_k**2)))*r_vect_hat_k 
            Rightside = (1- (3*(np.dot(r_vect_hat_j  ,r_vect_hat_k))/(r_dist_hat_j**2)))*r_vect_hat_j 
            
           
            
            Lindex = coef * (Leftside + Rightside)
            
            L_matrix.append(Lindex)
            
            
    L_matrix_nump = np.array(L_matrix).reshape((4,4,2))
    
    L_matrix_nump_X = L_matrix_nump[:,:,0]
    L_matrix_nump_Y = L_matrix_nump[:,:,1]
  
    
    
 
    
    N_Matrix_X = np.matmul(np.matmul(KI_matrix.transpose() , L_matrix_nump_X), KI_matrix)
    N_Matrix_Y = np.matmul(np.matmul(KI_matrix.transpose() , L_matrix_nump_Y), KI_matrix)
    
   
    I_matrix = np.array([ILeft,ITop,IRight,IBottom])
    
    Fx = ki * np.matmul(np.matmul(I_matrix, N_Matrix_X), I_matrix.transpose())
    Fy = ki * np.matmul(np.matmul(I_matrix, N_Matrix_Y), I_matrix.transpose())
    
   
  
   
    FXX.append(Fx*(1/pico))
    FYY.append(Fy*(1/pico)) # off by an order of magnetitude for  Fy according to paper
    


plt.figure()
#plt.quiver(arrow_start_X, arrow_start_Y ,FXX,FYY)#,scale=10, units = "xy", scale_units = "xy")
plt.quiver(measurement_location[0],measurement_location[1],FXX,FYY)
plt.title("Force Envelope at %s %s"%([measurement_location[0]*(1/micro),measurement_location[1]*(1/micro)],"um"))


plt.xlabel("FX")
plt.ylabel("FY")
plt.show()
##################################################################################################################
#%%







