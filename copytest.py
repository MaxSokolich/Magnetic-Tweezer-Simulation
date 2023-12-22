
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def Q_matrix(Nc, perm0, Ra, Ki, I):
    """
    args:
        Nc: number of turns on coil
        perm0: magnetic permeabilty of free space
        Ra: lumped reluctance from the pole top to the worspace center in the air
        Ki: the distribution matrix of the magnetic flux
        I: current matrix

    returns:
        Q: 4x1 matrix of charges [q1,q2,q3,q4] for a quadropole system  
    """
    Q_matrix = (Nc/(perm0*Ra)) * np.matmul(Ki , I)
    return Q_matrix



def reluctance(l_R, perm0, perm_R, A_R):
    """
    args:
       l_R: length of object
       perm0: magnetic permeability of free space
       perm_R: magnetic permeability of the object (mu = perm0*perm_R)
       A_R: cross section area of the object

    returns:
        Ra: magnetic relutance of an object. (similar to resistance)
    """
    R = l_R / (perm0 * perm_R * A_R)
    return R


def Ki_matrix():
    """
    args:
       None

    returns:
        Ki: magnetic curcuit distribution matrix of the magnetic flux
    """
    Ki = np.array([[3/4, -1/4, -1/4, -1/4],  
                      [-1/4, 3/4, -1/4, -1/4],
                      [-1/4, -1/4, 3/4, -1/4],
                      [-1/4, -1/4, -1/4, 3/4]])

    return Ki


def I_matrix():
    """
    args:
       None

    returns:
        I: current to each pole [I1,I2,I3,I4]
    """
    I = np.array([1,0,0,0])
    return I







def magnetic_charge(flux, perm0):
    """
    args:
        psi/flux: magnetic flux from the pole
        perm0: magnetic permeability of free space

    returns:
        q: the magnetic charge associated with that flux
    """
    q = flux/perm0
    return q


def magnetic_field(q, r, km, u):
    """
    args:
        q: magnetic point charge
        r: norm of the vector originating from the location of the magnetic charge to the magnetic particle
        km: constant defined by perm0/4pi
        u: unit directional vector of r;  u_ = r_ / r

    returns:
        B_: the magnetic field produced by the magnetic charges in the workspace. a 2D vector 
    """

    B = km * (q/r**2) * u
    return np.array(B)


def r(pole_location, measurement_location):
    """
    args:
        pole_location: 2D coordinates of pole location
        measurement_location: 2D coordinates of measurment location/ bot location

    returns:
        r_vector: the vector oroginating from the the location of the magnetic charge to the magnetic particle
        r_norm: the distance of r_vector
        u_vector: the unit normal vector of r_vector
    """
    x1 = pole_location[0]
    y1 = pole_location[1]
    
    x2 = measurement_location[0]
    y2 = measurement_location[1]
    
    r_vector = np.array([(x2-x1),(y2-y1)])
    r_norm = (np.sqrt((x2-x1)**2 + (y2-y1)**2))  
    u_vector = r_vector/r_norm

    return r_vector, r_norm, u_vector






def create_poles():
    pole1_location = np.array([0, 500*micro])
    pole2_location = np.array([500*micro,0])
    pole3_location = np.array([0,-500*micro])
    pole4_location = np.array([-500*micro,0])
    return [pole1_location,pole2_location,pole3_location,pole4_location] 


def find_field(pole_list, measurement_location):
    Bx = 0
    By = 0
    for j in range(len(pole_list)):

        r_vectorj, r_normj, u_vectorj = r(pole_list[j],measurement_location)
        Bj = magnetic_field(Q[j], r_normj, km, u_vectorj)
        
        Bx += Bj[0]
        By += Bj[1]

    return np.array([Bx,By])

def plot_magnetic_field():
    Pole_Pole_Distance = 1000 * micro
    X = np.linspace(-(Pole_Pole_Distance/2), (Pole_Pole_Distance/2), 50)
    Y = np.linspace(-(Pole_Pole_Distance/2),(Pole_Pole_Distance/2), 50)
    x,y = np.meshgrid(X, Y)

    fig, ax = plt.subplots()
    Bfield = find_field(pole_list,[x,y])
    
    #ax.quiver(x,y,Bfield[0],Bfield[1],color='r', scale = None)
    speed = np.sqrt((Bfield[0])**2+(Bfield[1])**2)
    ax.streamplot(x,y,Bfield[0],Bfield[1], color=speed,cmap='coolwarm',density = 2,norm = colors.LogNorm(vmin=speed.min(), vmax=speed.max() ))
    plt.show()


if __name__ == "__main__":
    milli = (10**(-3))
    micro = (10**(-6))
    nano = (10**(-9))
    pico = (10**(-12))
    #constants
    
    Nc = 400  #number of turns
    perm0 = ((4*np.pi)*(10**(-7)))
    km = (perm0/(4*np.pi))
    l_R = 0.405 * milli#length of air object
    perm_R = 1 #relative permeability of air
    A_R = 0.18 *milli *milli#cross section area of air object

    
    

    Ra = reluctance(l_R,perm0,perm_R, A_R) #check
    
    I = I_matrix()
    Ki = Ki_matrix()
    Q = Q_matrix(Nc,perm0,Ra, Ki, I) #check

    
    

    measurement_location = [0*micro,200*micro]
    
    pole_list = create_poles()

    plot_magnetic_field()
    
    

    
    
    
    
    
    
    
    
    
        
    
    
    



    

   
    

    
