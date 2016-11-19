import math
import numpy as np
import pylab as plt

def f(t,x,u):
    u1,u2 = u
    x1,x2,x3 = x
    f1 = -(3.8462*x1+0.2166*x3**2*math.sin(x2)-19.5258*math.cos(x2)*math.sin(x2)+0.00011352*x3*math.cos(x2)+0.1083*x1*x3*math.sin(x2)-0.1598*x1**2*math.cos(x2)**2*math.sin(x2)+0.1083*x1*x3*math.cos(x2)*math.sin(x2))/(0.1083*math.sin(x2)-0.3197*math.cos(x2)**2+1)
    f2 = x3                               
    f3 = -(0.000524*x3-90.1468*math.sin(x2)-9.7629*math.sin(x2)**2+5.6766*x1*math.cos(x2)+0.00005676*x3*math.sin(x2)-0.0799*x1**2*math.cos(x2)*math.sin(x2)**2-0.7379*x1**2*math.cos(x2)*math.sin(x2)+0.3197*x3**2*math.cos(x2)*math.sin(x2)+0.0160*x1*x3*math.cos(x2)*math.sin(x2)+0.1598*x1*x3*math.cos(x2)**2*math.sin(x2))/(0.1083*math.sin(x2)-0.3197*math.cos(x2)**2+1)

    g1 = 76.9231/(0.1083*math.sin(x2)-0.3197*math.cos(x2)**2+1)
    g2 = 0
    g3 = 113.5308*math.cos(x2)/(0.1083*math.sin(x2)-0.3197*math.cos(x2)**2+1)
    
    F1 = f1+g1*u1
    F2 = f2+g2*u1
    F3 = f3+g3*u1
    return np.array([F1,F2,F3])


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
        x1,x2,x3 = x_[-1]+dt*f(j*dt,x_[-1],u_[-1])
	'''
        x2 = x2 % (2*math.pi)
        if x2 > math.pi:
            x2 -= 2*math.pi
	'''
        x_.append(np.array([x1,x2,x3]))
        j += 1
    return np.array(t_),np.array(x_)



t1 = 10
xs = np.array([0,math.pi/4.0,0])
us = [0]
t_,x_ = sim(f,t1,xs,us,dt=1e-4)

lineObjects = plt.plot(t_,x_,'-')
plt.xlabel('time (sec)')
plt.ylabel('state (init =[0,math.pi/4.0,0])')
plt.legend(iter(lineObjects),('x1 = dAlpha','x2 = Beta','x3 = dBeta'))
plt.show()
