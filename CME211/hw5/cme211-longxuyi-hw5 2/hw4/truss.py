import sys
import os
import math
import numpy as np
import warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
import scipy.sparse.linalg
from scipy.sparse.linalg import MatrixRankWarning


class Truss:

    def __init__ (self, joints_file, beams_file, plot_file):       
        self.beamforce(joints_file, beams_file, plot_file)
 

    def beamforce(self, joints_file, beams_file, plot_file):
        '''this function calculate force on each beam
        and raise error if the system is either over or under determined'''

        x_pos=[]
        y_pos=[]
        self.joints=[]
        self.jointId = []
        self.joints_dict = {}
        
        f = open(joints_file, "r")
        lines = f.readlines()[1:]
        for line in lines:
                id,x,y,Fx,Fy,zerodisp = int(line.split()[0]),float(line.split()[1]), float(line.split()[2]), float(line.split()[3]), float(line.split()[4]), float(line.split()[5])
                x_pos.append(x)
                y_pos.append(y)
                self.joints.append([x,y])
                self.joints_dict[id] = [x,y,Fx, Fy, zerodisp]
        f.close()

        joint_a=[]
        joint_b=[]
        self.beams=[]
        self.beams_dict = {}
          
        f = open(beams_file, "r")
        lines = f.readlines()[1:]
        for line in lines:
                id, ja,jb = int(line.split()[0]),int(line.split()[1]), int(line.split()[2])
                joint_a.append(x)
                joint_b.append(y)
                self.beams.append([self.joints[ja-1], self.joints[jb-1]])  
                self.beams_dict[id] = [ja, jb]
                #self.beamsId.append(id)    
        f.close()

        #make a optional plot of the geometry
        if plot_file != "":
            self.PlotGeometry(plot_file)

        #check if number of unknowns matches number of equations
        n = 0
        for key in self.joints_dict:
            if self.joints_dict[key][4]==1:
                n= n+2
        n = n +len(self.beams_dict)

        if n != 2*len(self.joints_dict):
            raise RuntimeError("Truss geometry not suitable for static equilibrium analysis")

        #construct a COO cefficient matrix
        row = []
        col = []
        data =[]
        rigid_pt =0

        #loop through each joint
        for i in range(2*len(self.joints_dict)):          
            #loop thru x components
            if i%2 ==0:
                #loop through columns
                for j in range(n):
                    
                    #for beam colomns
                    if j <= len(self.beams_dict)-1:

                        #if the beam links to the joint 
                        if (i/2+1 in self.beams_dict[j+1]):
                            row.append(i)
                            col.append(j)
                            
                            joint_a = self.beams_dict[j+1][0]
                            joint_b = self.beams_dict[j+1][1]

                            ax=self.joints_dict[joint_a][0]
                            ay=self.joints_dict[joint_a][1]
                            bx=self.joints_dict[joint_b][0]
                            by=self.joints_dict[joint_b][1]

                            if int(i/2+1) == joint_a: 
                                val = (ax -bx)/((bx-ax)**2+(by-ay)**2)**0.5
                            else:
                                val = (bx -ax)/((bx-ax)**2+(by-ay)**2)**0.5
                            
                            data.append(val)
                    
                #values for rigid points
                if self.joints_dict[i/2+1][4]==1:

                                rigid_pt += 1
                                row.append(i)
                                col.append(len(self.beams_dict)-1+rigid_pt*2-1)
                                data.append(1)

            #loop through y components
            else:
                for j in range(n):   
                    #for beam colomns
                    if j <= len(self.beams_dict)-1:

                        #if the beam links to the joint 
                        if ((i-1)/2+1 in self.beams_dict[j+1]):
                            row.append(i)
                            col.append(j)
                            
                            joint_a = self.beams_dict[j+1][0]
                            joint_b = self.beams_dict[j+1][1]

                            ax=self.joints_dict[joint_a][0]
                            ay=self.joints_dict[joint_a][1]
                            bx=self.joints_dict[joint_b][0]
                            by=self.joints_dict[joint_b][1]

                            if int((i-1)/2+1) == joint_a: 
                                val = (ay -by)/((bx-ax)**2+(by-ay)**2)**0.5
                            else:
                                val = (by -ay)/((bx-ax)**2+(by-ay)**2)**0.5
                            
                            data.append(val)
                                        #values for rigid points
                    
                if self.joints_dict[(i-1)/2+1][4]==1:
                                row.append(i)
                                col.append(len(self.beams_dict)-1+rigid_pt*2)
                                data.append(1)
        
        #construct an array to store external forces
        ext_force =np.zeros((n,1))

        # create external forces matrix
        for i in self.joints_dict:
            if self.joints_dict[i][2] != 0:
                ext_force[2*(i-1),0] = self.joints_dict[i][2]*(-1)
            if self.joints_dict[i][3] != 0:
                ext_force[2*(i-1)+1,0] = self.joints_dict[i][3]*(-1)
        b = scipy.sparse.csr_matrix(ext_force)
        a = coo_matrix((data, (row, col)), shape=(n, n)).tocsr()
        self.x = scipy.sparse.linalg.spsolve(a,b)

        #raise error if the matrix is unsolvable
        try:
            warnings.filterwarnings('error')
            self.x = scipy.sparse.linalg.spsolve(a,b)
        except MatrixRankWarning:
            raise RuntimeError("The linear system is not solvable, unstable truss?")


    def PlotGeometry(self, beam_plot):
        '''this function creates a geometry of beam structure'''
        for i in range(len(self.beams)):
            x_values = [self.beams[i][0][0], self.beams[i][1][0]]
            y_values = [self.beams[i][0][1], self.beams[i][1][1]]
            plt.plot(x_values, y_values, 'bo', linestyle='solid')   

        plt.gca().margins(0.05,0.05)
        plt.axes().set_aspect('equal', 'datalim')
        plt.savefig('{}'.format(beam_plot))
   
        
    def __repr__(self):
        '''this string representation function prints computed beam forces'''
        header1 = "Beam   Force\n"
        header2 = '------------\n'
        output = header1 + header2
        for i in range(len(self.beams_dict)):
            string = "{}     {:.3f}\n".format(i+1, self.x[i])
            output += string
        return output          