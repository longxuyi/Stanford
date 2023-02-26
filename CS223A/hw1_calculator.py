import numpy as np
import math
 
np.set_printoptions(precision=3)
# input two matrices
# mat1 = ([0, 0, 1],[-0.5 ,0.866, 0],[0.866, 0.5, 0])
# mat2 = ([0.866, -0.5, 0],[0.5, 0.866, 0],[0,0,1])

# cp = math.cos(math.pi/3)
# sp = math.sin(math.pi/3)
# ct = math.cos(math.pi/4)
# st = math.sin(math.pi/4)

# mat1 = ([cp, 0, sp],[0 ,1, 0],[-sp, 0, cp])
# mat2 = ([1, 0, 0],[0, ct, -st],[0,st,ct])
 
# # This will return dot product
# rmat = np.dot(mat1,mat2)

# vec = ([3],[1],[2])

# rvec = np.dot(rmat, vec)

# # print resulted matrix
# print(rvec)

# problem 3 (b)
cp = math.cos(math.pi/6)
sp = math.sin(math.pi/6)
ct = math.cos(math.pi/4)
st = math.sin(math.pi/4)


mat1 = ([1, 0, 0],[0, ct, -st],[0,st,ct])
mat2 = ([cp, 0, sp],[0 ,1, 0],[-sp, 0, cp])
mat3 = ([ct, -st, 0],[st, ct, 0],[0,0,1])

res = np.dot(mat3,mat2)

rmat = np.dot(res, mat1)

print("3b:")
print(rmat)

# problem 3c
cp = math.cos(math.pi/2)
sp = math.sin(math.pi/2)
ct = math.cos(math.pi/4)
st = math.sin(math.pi/4)


mat1 = ([1, 0, 0],[0, ct, -st],[0,st,ct])
mat2 = ([cp, 0, sp],[0 ,1, 0],[-sp, 0, cp])
mat3 = ([ct, -st, 0],[st, ct, 0],[0,0,1])

res = np.dot(mat3,mat2)

rmat = np.dot(res, mat1)

# print(mat1)
# print(mat2)
# print(mat3)
print("3c:")
print(rmat)


# PROBLEM 4b
rmat = np.array([[0.5, -0.866, 0],
              [0.866, 0.5, 0],
              [0, 0, 1]])
  
# Calculating the inverse of the matrix
inv = np.linalg.inv(rmat)

print(inv)

vec0 = ([2],[1],[2])

#inverse of translation vector
invec = -1*np.dot(inv, vec0)

print(invec)



inv = np.array([[0.5, 0.866, 0, -1.866],
              [-0.866, 0.5, 0, 1.232],
              [0, 0, 1, -2],
              [0, 0, 0, 1]])


vec = ([3],[1],[2],[1])
  
res = np.dot(inv,vec)


print(res)

#problem 4c

# #beta
# print(math.atan2(0, 1))

# #alpha
# print(math.degrees(math.atan2(0.866, 0.5)))

# #gemma
# print(math.atan2(0, 1))

# #angle-axis representation

# print(math.degrees(math.acos( (0.5+0.5 +1 -1)/2)))


r = np.array([[0.866, 0.25, -0.433],
              [0, 0.866, 0.5],
              [0.5, -0.433, 0.75]])


svas = ([-0.25],[1.5],[-0.1])

spa = ([1.5],[4],[0])

gvsg = ([-1],[5],[2])

ss = ([0],[0],[0.35])

gs = np.array(np.dot(r,ss))

mat = np.array(np.dot(r,spa))

print("mat:")
print(mat)

print("gs:")
print(gs)

result = gvsg + np.dot(r,svas) + np.cross(gs.flatten(), mat.flatten()).reshape(3,1)

print("result:")
print(result)
