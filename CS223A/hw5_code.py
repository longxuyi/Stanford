
import matplotlib.pyplot as plt
import numpy as np
import math

#define start, via, and end points
start = [2, 2]
via = [6, 5]
end = [7.5, 10] 

#define speed at start, via, and end points
velo_start =0
velo_via = 1
velo_end = 0

#define time at via point
tv = 3

#calculate a0,a1,a2,a3 for each cubic polynomial 
def solver(pos_0, pos_f, velo_0, velo_f, tf ):
    p0 = pos_0
    p1 = velo_0
    p2 = 3/pow(tf,2)*(pos_f - pos_0) - 2/tf * velo_0 -1/tf * velo_f
    p3 = -2/pow(tf,3)*(pos_f - pos_0) + 1/pow(tf,2)*(velo_0 + velo_f)

    # print("u(t) = {} + {} * t + {} * t^2 + {} * t^3".format(round(p0,3),round(p1,3),round(p2,3),round(p3,3)))
    return p0, p1, p2, p3

#unknowns for x
a0,a1,a2,a3 = solver(start[0], via[0], velo_start, velo_via, tv)
b0,b1,b2,b3 = solver(via[0], end[0], velo_via, velo_end, 6-tv)

#unknowns for y
c0,c1,c2,c3 = solver(start[1], via[1], velo_start, velo_via, tv)
d0,d1,d2,d3 = solver(via[1], end[1], velo_via, velo_end, 6-tv)

#piece wise function for x
def u(t):
    if t < tv:
        return a0 + a1 * t + a2 * pow(t,2) + a3 * pow(t,3)
    if t >= tv:
        return b0 + b1 * (t-tv) + b2 * pow(t-tv,2) + b3 * pow(t-tv,3)

#piece wise function for y
def v(t):
    if t < tv:
        return c0 + c1 * t + c2*math.pow(t,2) + c3 * math.pow(t,3)
    if t >= tv:
        return d0 + d1 * (t-tv) + d2 * math.pow((t-tv),2) + d3 * math.pow((t-tv),3)

t1 = np.linspace(0,tv,100*tv)
t2 = np.linspace(tv,6,100*(6 -tv))

x1= []
x2 =[]
y1= []
y2 =[]

for i in range(len(t1)):
    x1.append(u(t1[i]))
    y1.append(v(t1[i]))

for i in range(len(t2)):
    x2.append(u(t2[i]))
    y2.append(v(t2[i]))

#plot 2d path
x = x1 + x2
y = y1 + y2
t = np.linspace(0,6,600)

slope_start = (y[1] -y[0])/(x[1] -x[0])
slope_via = (y[301] -y[299])/(x[301] -x[299])
slope_end = (y[599] -y[598])/(x[599] -x[598])

print("start: dy/dt = {}".format(slope_start))
print("via: dy/dt = {}".format(slope_via))
print("end: dy/dt = {}".format(slope_end))

plt.plot(t,x, label='x')
plt.plot(t,y, label='y')
plt.xlabel("t")
plt.ylabel("x and y")
plt.title("x and y positions vs time")
plt.legend()
plt.grid()
plt.show()


plt.plot(x, y, "r", linewidth=2)

#plot circles
c1 = plt.Circle((6.5, 1.5),2.5 ,fill = False )
c2 = plt.Circle((4, 7),2 ,fill = False )
plt.gca().add_artist(c1)
plt.gca().add_artist(c2)

#plot start, via, and end point
plt.plot(start[0], start[1], marker="o", markersize=10, markeredgecolor="red", markerfacecolor="green")
plt.plot(via[0], via[1], marker="+", markersize=10, markeredgecolor="black", markerfacecolor="black")
plt.plot(end[0], end[1], marker="*", markersize=15, markeredgecolor="red", markerfacecolor="purple")

plt.axis('equal')
plt.xticks(np.arange(0, 12, 1))
plt.yticks(np.arange(-3, 15, 1))
plt.xlabel("X(t)")
plt.ylabel("Y(t)")
plt.title("Motion Planning")

plt.grid()
plt.show()
