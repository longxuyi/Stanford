import sys
import glob
import math
import os

class Airfoil:
    '''This class handles the loading and processing of data 
    associated with an airfoil at multiple angles of attack'''
    
    def __init__ (self, naca_folder):
        '''this function process all files inside the directory input and
        summarize test data'''
        
        geometry = self.load_geometry(naca_folder)
        alpha = self.load_alpha(naca_folder)
        self.find_lift_coeff(geometry, alpha)
        self.find_stag_point(geometry, alpha)
        self.generate_summary()

    
    def load_geometry(self, naca_folder):
        '''create a list containing (x,y) coordiniates for the panel'''

        self.geometry = []
        self.path = naca_folder

        isExist = os.path.exists(self.path)
        if not isExist:
            raise RuntimeError("not a valid directory")

        try:
            f = open(self.path + "/xy.dat", 'r')
        except IOError:
            raise RuntimeError("data file cannot be found")

        try:
            #skip the headerline
            lines = f.readlines()[1:]
        except:
            raise RuntimeError("errors in the data file")
  
        for line in lines:
                x,y = float(line.split()[0]), float(line.split()[1])
                self.geometry.append([x,y])
        f.close()
   
        return self.geometry


    def load_alpha(self, naca_folder):
        """create a dictionary to store Cp values for each angle. 
        With each angle as a key and corresponding a list of Cp 
        """
        self.Cp = []
        self.alpha = {}
        dir_list = os.listdir(naca_folder)

        #find alpha files and create a dicctionary with alpha 
        #as a key and list of Cp as value 
        for file in dir_list:
            #find alpha file
            if file.startswith('alpha'):
                key = float(file[5:-4])
                #initialize a list of Cp value
                Cp =[]
                f = open(self.path + "/" + file, 'r')
                lines = f.readlines()[1:]

                for line in lines:
                    Cp_value = float(line.split()[0])
                    Cp.append(Cp_value)
                
                f.close()
            
                self.alpha[key] = Cp 

        return self.alpha


    def find_lift_coeff(self, geometry, alpha):

        '''calculate lift coefficient'''

        self.lift_coefficient = []
        
        #cord length calculation
        x_value = []

        for i in range(len(geometry)):

            x_value.append(geometry[i][0])

        x_max = max(x_value)
        max_index = x_value.index(x_max)
        x_trailing = x_max
        y_trailing = geometry[max_index][1]
        chord_length = (x_trailing**2 +y_trailing**2)**0.5

        for angle in alpha:
            cos_alpha = math.cos(math.pi/2*float(angle)/90)
            sin_alpha = math.sin(math.pi/2*float(angle)/90)
            
            cx = 0
            cy = 0

            # calculate lift coefficient for each alpha
            for i in range(len(geometry)-1):

                delta_cx = -1*((geometry[i+1][1] - geometry[i][1])* alpha[angle][i])/chord_length
                delta_cy = ((geometry[i+1][0] - geometry[i][0])* alpha[angle][i])/chord_length

                cx += delta_cx
                cy += delta_cy

            lc = round((cy*cos_alpha - cx*sin_alpha),4)

            self.lift_coefficient.append(lc)


    def find_stag_point(self, geometry, alpha):
        '''This function calculates the stagnation points 
        at different alpha angles and return a list of
        stagnation points'''

        self.stag_points =[]

        self.stag_Cp = []

        for angle in alpha:

            max_Cp = max(alpha[angle])

            max_index = alpha[angle].index(max_Cp)

            #find average of two points defining the pannel at the highest Cp
            stag_x = round((geometry[max_index+1][0]+ geometry[max_index][0])/2,4)
            stag_y = round((geometry[max_index+1][1]+ geometry[max_index][1])/2,4)

            stag = [stag_x, stag_y]

            self.stag_points.append(stag)

            self.stag_Cp.append(round(max_Cp,4))
      
        
    def generate_summary(self):
        ''' Generate a summary of lift coefficient 
        and stagnation points at different alpha angles'''

        list_summary = []

        #convert key set to list 
        alpha_angles = sorted(list(self.alpha.keys()))

        for i in range(len(alpha_angles)):

            list_summary.append([alpha_angles[i], self.lift_coefficient[i],self.stag_points[i], self.stag_Cp[i]])
        
        self.list_summary = sorted(list_summary)


    def __repr__(self):
        '''this string representation function prevents main.py to 
        print Airfoil as a string rather than an object.'''

        header1 = "alpha       cl          stagnation point\n"
        header2 = '------    ------     ----------------------\n'
        output = header1 + header2

        for self.summary in self.list_summary:
            string = "{}     {}     ({},{})      {}\n".format(self.summary[0], self.summary[1],
            self.summary[2][0], self.summary[2][1], self.summary[3])
            output += string
        
        return output
