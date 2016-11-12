import math
import numpy as np
import pylab as plt

k = 0.5

def f(t,x,u):
    x1,x2,v1,v2,theta = x
    u1,u2 = u
    return np.array([v1,v2,u1*math.cos(theta)-k*v1,u1*math.sin(theta)-k*v2,u2])

def h(t,x,u):
    x1,x2,v1,v2,theta = x
    u1,u2 = u
    return np.array([x1,x2,theta])

def u1(t):
    return np.array([0,0])

def sim(f, t, x, u, dt=1e-4):
    j, t_, x_, u_ = 0, [0], [x], [u]
    while j*dt < t:
        t_.append((j+1)*dt)
        u_.append(u1(j*dt))
        x_.append(x_[-1]+dt*f(j*dt,x_[-1],u_[-1]))
        j += 1
    return np.array(t_),np.array(x_)

#x0 = [1,1,0,0,2]
#u0 = [0,0]
#print 'f(x0) = ', f(0.,x0,u0)


t1 = 10
xs = [2,1,-1,1,-1.57]
us = [0,0]
t_,x_ = sim(f,t1,xs,us,dt=4e-2)

lineObjects = plt.plot(t_,x_,'-')
plt.xlabel('time (sec)')
plt.ylabel('state (init =[2,1,-1,1,-1.57])')
plt.legend(iter(lineObjects),('x1','x2','v1','v2','theta'))
plt.show()










