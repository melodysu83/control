import numpy as np
import pylab as plt
import math

k = 1

def f(t,x,u):
    x1,x2,v1,v2,theta = x
    u1,u2 = u
    return np.array([v1,v2,u1*math.cos(theta)-k*v1,u1*math.sin(theta)-k*v2,u2])

def h(t,x,u):
    x1,x2,v1,v2,theta = x
    u1,u2 = u
    return np.array([x1,x2,theta])

def get_u(t):
    u1 = 1
    u2 = 1/(1+0.1*t)
    return np.array([u1,u2])

def sim(f, t, x, u, dt=1e-4):
    j, t_, x_, x1_, x2_, u_ = 0, [0], [x], [x[0]], [x[1]], [u]
    while j*dt < t:
        t_.append((j+1)*dt)
        u_.append(get_u(j*dt))
        x1,x2,v1,v2,theta = x_[-1]+dt*f(j*dt,x_[-1],u_[-1])
        theta = theta % (2*math.pi)
        if theta > math.pi:
            theta -= 2*math.pi
        x1_.append(x1)
        x2_.append(x2)
        x_.append(np.array([x1,x2,v1,v2,theta]))
        j += 1
    return np.array(t_),np.array(x_),np.array(x1_),np.array(x2_)

total_t = 60
x0 = [1,1,-1,-1,2]
u0 = [0,0]
print 'f(x0) = ', f(0.,x0,u0)

t_,x_,x1_,x2_ = sim(f,total_t,x0,u0,dt=4e-2)

plt.figure()
lineObjects = plt.plot(t_,x_,'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.legend(iter(lineObjects),('x1','x2','v1','v2','theta'))
plt.title('State-Time plot')

plt.figure()
plt.plot(x1_,x2_,'-')
plt.xlabel('x1')
plt.ylabel('x2')
plt.title('x1-x2 plot')

plt.figure()
for i in range(5):
    xinit = x0+np.array([0,0,i,0,10*i*math.pi/180.0])
    t_,x_,x1_,x2_ = sim(f,total_t,xinit,u0,dt=4e-2)
    plt.plot(t_,x1_,'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.title('x1-time plot (with different initial condition)')

plt.figure()
for i in range(5):
    xinit = x0+np.array([0,0,i,0,10*i*math.pi/180.0])
    t_,x_,x1_,x2_ = sim(f,total_t,xinit,u0,dt=4e-2)
    plt.plot(t_,x2_,'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.title('x2-time plot (with different initial condition)')

plt.figure()
for i in range(5):
    xinit = x0+np.array([0,0,i,0,10*i*math.pi/180.0])
    t_,x_,x1_,x2_ = sim(f,total_t,xinit,u0,dt=4e-2)
    plt.plot(t_,x_[:,2],'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.title('x3-time plot (with different initial condition)')

plt.figure()
for i in range(5):
    xinit = x0+np.array([0,0,i,0,10*i*math.pi/180.0])
    t_,x_,x1_,x2_ = sim(f,total_t,xinit,u0,dt=4e-2)
    plt.plot(t_,x_[:,3],'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.title('x4-time plot (with different initial condition)')

plt.figure()
for i in range(5):
    xinit = x0+np.array([0,0,i,0,10*i*math.pi/180.0])
    t_,x_,x1_,x2_ = sim(f,total_t,xinit,u0,dt=4e-2)
    plt.plot(t_,x_[:,4],'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.title('x5-time plot (with different initial condition)')

plt.show()
