# x -> infected
# y -> susceptible 
# z -> removed(both died and recovered patients)
# x(t)+ y(t)+z(t)=N where N is population of the country 

# x(0)=x_0
# y(0)=y_0 = N - x_0
# z(0)=0

# dx/dt = Axy - Bx
# dx/dt = -Axy 
# dz/dt = Bx
# A -> infection rate 
# B -> removal rate 
# p = B/A 
# According the theorem of the Epidemiology-> Epidemic will break if in the population if y_0 > p

import numpy as np
import matplotlib.pyplot as mp
import math as m
# Euler's formula -> x[t]= x[t-1]+ f'(x)del_t
def euler_integrate_2(a,b,x0,y0, del_t,A,B):
    n = int((b-a)/del_t)-1
    x = np.zeros(n+1)
    y = np.zeros(n+1)
    x[0]= x0
    y[0]=y0
    for i in range(1, n+1):
        x[i]= x[i-1]+(A*x[i-1]*y[i-1]-B*x[i-1]) *del_t
        y[i]= y[i-1]-(A*x[i-1]*y[i-1]) *del_t
    return y,x    

def euler_integrate_3(a, b,x0,y0,z0, del_t,A,B):

    n = int((b-a)/del_t) -1
    x = np.zeros(n+1)
    y = np.zeros(n+1)
    z = np.zeros(n+1)
    x[0] = x0
    y[0]= y0
    z[0]=z0

    for i in range(1, n+1):
        x[i] = x[i-1] + del_t*(A*x[i-1]*y[i-1] - B*x[i-1])
        y[i] = y[i-1] - del_t*A*x[i-1]*y[i-1]
        z[i] = z[i-1] + del_t*B*x[i-1]

    return x,y,z

#Part 1 : plottig y/x for different parameters
# disease will break 
N = 1e5
x_0 = 10
A =10
B=1
y_0= N -x_0
z_0 =0
initVal =0
del_t = 0.0000005 # t is in time 
y_euler, x_euler = euler_integrate_2(initVal,100000*del_t, x_0, y_0,del_t,A,B)
mp.figure()
mp.plot(y_euler,x_euler,linewidth=1.5)
mp.grid(which='major', color='#DDDDDD', linewidth = 1)
mp.xlabel('y(t)')
mp.ylabel('x(t)')
mp.title('Plot of x(t)vs y(t)')
mp.legend(["B/A=10"])


# In this graph the disease will not break 
N = 1e5
x_0 = 10
A =0.0000001
B=1
y_0= N -x_0
z_0 =0
initVal =0
del_t = 0.0005 # t is in time 
y_euler, x_euler = euler_integrate_2(initVal,1000*del_t, x_0, y_0,del_t,A,B)
mp.figure()
mp.plot(y_euler,x_euler,linewidth=1.5)
mp.grid(which='major', color='#DDDDDD', linewidth = 1)
mp.xlabel('y(t)')
mp.ylabel('x(t)')
mp.title('Plot of x(t)vs y(t)')
mp.legend(["B/A>y_0"])



#Part 2 : plottig x,y,z vs t
# disease will break 
N = 100000
x_0 = 10
A =10
B=1
y_0= N -x_0
z_0 =0
initVal =0
del_t = 0.00000005 # t is in time 
x_euler, y_euler,z_euler = euler_integrate_3(initVal,1000*del_t, x_0, y_0,z_0,del_t,A,B)

timeRange = np.arange(initVal,1000*del_t,del_t)
mp.figure()
mp.plot(timeRange,x_euler,linewidth=1.5)
mp.plot(timeRange,y_euler,linewidth=1.5)
mp.plot(timeRange,z_euler,linewidth=1.5)
mp.xlabel('t')
mp.ylabel('x(t), y(t)  and z(t)')
mp.legend(["x(t)", "y(t)" ,"z(t)"])
mp.title('Plot of x(t),y(t) & z(t)')
mp.grid()

mp.figure()
mp.plot(timeRange,x_euler,linewidth=1.5)
mp.xlabel('t')
mp.ylabel('x(t)')
mp.legend(["x(t)"])
mp.title('Plot of x(t)')
mp.grid()

mp.figure()
mp.plot(timeRange,y_euler,linewidth=1.5, color='orange')
mp.xlabel('t')
mp.ylabel('y(t)')
mp.legend([ "y(t)" ])
mp.title('Plot of y(t)')
mp.grid()

mp.figure()
mp.plot(timeRange,z_euler,linewidth=1.5,color='green')
mp.xlabel('t')
mp.ylabel('z(t)')
mp.legend(["z(t)"])
mp.title('Plot ofz(t)')
mp.grid()

x_euler_log = np.log(x_euler)
z_euler_log =np.log(z_euler)
mp.figure()
mp.plot(timeRange,x_euler_log,linewidth=1.5)
mp.plot(timeRange,z_euler_log,linewidth=1.5)
mp.xlabel('t')
mp.ylabel('x(t) z(t)')
mp.legend(["ln x(t)" ,"ln z(t)"])
mp.title('Plot of ln x(t)& ln z(t)')
mp.grid()

#Another parameters
N = 100000
x_0 = 100
A =0.000001
B=1
y_0= N -x_0
z_0 =0
initVal =0
del_t = 0.0005 # t is in time 
x_euler, y_euler,z_euler = euler_integrate_3(initVal,1000000*del_t, x_0, y_0,z_0,del_t,A,B)

timeRange = np.arange(initVal,1000000*del_t,del_t)
mp.figure()
mp.plot(timeRange,x_euler,linewidth=1.5)
mp.plot(timeRange,y_euler,linewidth=1.5)
mp.plot(timeRange,z_euler,linewidth=1.5)
mp.xlabel('t')
mp.ylabel('x(t), y(t)  and z(t)')
mp.legend(["x(t)", "y(t)" ,"z(t)"])
mp.title('Plot of x(t),y(t) & z(t)')
mp.grid()

mp.show()
