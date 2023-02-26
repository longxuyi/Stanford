import numpy as np
import math
 
#problem 1b DH parameter table
l2 = 0.1
d3 =0.25
theta_1 =0
theta_2 = math.pi/2

DH = [[0, 0, 0, theta_1],
      [0, math.pi/2, 0,theta_2],
      [l2, -math.pi/2,d3, 0]]

num_rows = np.shape(DH)[0]


print(num_rows)

#initialize T
T = np.array([[1,0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]])

for i in range(num_rows):
    np.set_printoptions(precision=3)

    #cos and sin alpha (i-1)
    ca = math.cos(DH[i][1])
    sa = math.sin(DH[i][1])
    
    #alpha (i-1)
    a = DH[i][0]

    #cos and sin theta i
    ct = math.cos(DH[i][3])
    st = math.sin(DH[i][3])

    #di
    d = DH[i][2]

    mat1 = np.array([[1,0, 0, 0],
            [0, ca, -sa, 0],
            [0,sa, ca, 0],
            [0, 0, 0, 1]])
    
    #print(np.around(mat1,3))

    mat2 = np.array([[1,0, 0, a],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]])
    
    #print(np.around(mat2,3))

    mat3 = np.array([[ct,-st, 0, 0],
            [st, ct, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]])

    #print(np.around(mat3,3))

    mat4 = np.array([[1,0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, d],
            [0, 0, 0, 1]])

    #print(np.around(mat4,3))
    
    temp =  np.dot(mat1, np.dot(mat2, np.dot(mat3,mat4)))

    print(np.around(temp,3))

    T = np.dot(T, temp)


print("transform matrix from 3 to 0:")

vec3 = ([0],[0],[0.1], [1])

vec0 = np.dot(T, vec3)

print(np.around(T,3))

# print("vector p in frame 0:")
# print(np.around(vec0,3))
    

    
# #problem 2b DH parameter table
# l1 = 1
# l3 =1
# d2 = 2
# theta_1 = math.pi/6
# theta_3 = -math.pi/3
# theta_4 =  math.pi/3


# DH = [[0, 0, 0, theta_1],
#       [l1, 0, d2,math.pi/2],
#       [0, math.pi/2, 0, theta_3],
#       [l3, math.pi/2, 0,theta_4]]

# num_rows = np.shape(DH)[0]

# #initialize T
# T = np.array([[1,0, 0, 0],
#             [0, 1, 0, 0],
#             [0, 0, 1, 0],
#             [0, 0, 0, 1]])

# for i in range(num_rows):
#     print(i)
#     np.set_printoptions(precision=3)

#     #cos and sin alpha (i-1)
#     ca = math.cos(DH[i][1])
#     sa = math.sin(DH[i][1])
    
#     #alpha (i-1)
#     a = DH[i][0]

#     #cos and sin theta i
#     ct = math.cos(DH[i][3])
#     st = math.sin(DH[i][3])

#     #di
#     d = DH[i][2]

#     mat1 = np.array([[1,0, 0, 0],
#             [0, ca, -sa, 0],
#             [0,sa, ca, 0],
#             [0, 0, 0, 1]])
#     #print(np.around(mat1,3))

#     mat2 = np.array([[1,0, 0, a],
#             [0, 1, 0, 0],
#             [0, 0, 1, 0],
#             [0, 0, 0, 1]])
#     #print(np.around(mat2,3))

#     mat3 = np.array([[ct,-st, 0, 0],
#             [st, ct, 0, 0],
#             [0, 0, 1, 0],
#             [0, 0, 0, 1]])
#     #print(np.around(mat3,3))

#     mat4 = np.array([[1,0, 0, 0],
#             [0, 1, 0, 0],
#             [0, 0, 1, d],
#             [0, 0, 0, 1]])
#     #print(np.around(mat4,3))
    
#     #temp =  np.dot(mat1, np.dot(mat2, np.dot(mat3,mat4)))

#     temp =  np.dot(np.dot(np.dot(mat1,mat2), mat3),mat4 )

#     print(np.around(temp,3))

#     T = np.dot(T, temp)

#     print(np.around(T,3))

# #print(np.around(T,3))

# #PROBLEM 2c
# T30 = np.array([[-0.25,-.433, .866, 0],
#             [.433, .75, .5, .5],
#             [-.866,.5, 0, 2],
#             [0, 0, 0, 1]])

# T34 = np.array([[ 0.5,  -0.866,  0,     1   ],
#                   [ 0,     0,    -1,     0   ],
#                   [ 0.866,  0.5,    0,     0   ],
#                   [ 0,    0,     0,     1   ]])

# T40 = np.dot(T30, T34)

#print(np.around(T40,3))