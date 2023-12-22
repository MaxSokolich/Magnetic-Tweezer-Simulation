import numpy as np
import matplotlib.pyplot as plt




class Pole:
    def __init__(self, pole_pos, measurement_location):
        self.pole_pos = pole_pos
        self.measurement_location = measurement_location

      
    def r_vector(self):
        '''
        returns a 2D array vector in direction of microbot from desired pole
        '''
        x1 = self.pole_pos[0]
        y1 = self.pole_pos[1]
        
        x2 = self.measurement_location[0]
        y2 = self.measurement_location[1]
        r_vector = np.array([(x2-x1),(y2-y1)])
    
        return r_vector
    

    def r_distance(self):
        '''
        returns a scalar distance desired pole to microbot
        not the normalized sistance r/L
        '''
        x1 = self.pole_pos[0]
        y1 = self.pole_pos[1]
        
        x2 = self.measurement_location[0]
        y2 = self.measurement_location[1]
        
        r_dist = (np.sqrt((x2-x1)**2 + (y2-y1)**2))

        return r_dist
    

    def create_pole(self, ax):
        pole_dist = self.r_distance()
        pole_vector = self.r_vector()/pole_dist

        ax.quiver(self.pole_pos[0],self.pole_pos[1], pole_vector[0], pole_vector[1], scale = None)









    


    


        










if __name__ == "__main__":

    fig, ax = plt.subplots()


    measurement_location = [100,200]







    Pole1 = Pole([0,400], measurement_location).create_pole(ax)
    Pole2 = Pole([400,0], measurement_location).create_pole(ax)
    Pole3 = Pole([0,-400], measurement_location).create_pole(ax)
    Pole4 = Pole([-400,0], measurement_location).create_pole(ax)

   


    ax.scatter(measurement_location[0],measurement_location[1],s=100,facecolors = "k")

    plt.show()
    

