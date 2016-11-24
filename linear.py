import numpy as np
import pylab as plt
import math
import copy

K = 1

A = np.matrix([[0,0,1,0,0],[0,0,0,1,0],[0,0,-K,0,0],[0,0,0,-K,0],[0,0,0,0,0]])
B = np.matrix([[0,0],[0,0],[1.0/math.sqrt(2),0],[1.0/math.sqrt(2),0],[0,1]])

total_t = 15 

ueq = np.array([0,0])
uo  = np.array([0,0])

Xeq  = np.array([1,1,0,0,math.pi/4.0])
Xo_1 = np.array([1,2,-0.5,0.5,1]) 
Xo_2 = np.array([3,4,0.5,-0.5,2]) 


def addM(a,b):
    res = []
    for i in range(len(a)):
        row = []
        for j in range(len(a[0])):
            row.append(a[i][j]+b[i][j])
        res.append(row)
    return np.array(res)   

def get_u(t, index):
    
    if index == 1:
        u1 = math.sin(t)
        u2 = math.cos(t)
        return np.array([u1,u2])
    elif index == 2:
        u1 = math.cos(t)
        u2 = math.sin(t)
        return np.array([u1,u2])
    else:
        return np.array([math.cos(t)+math.sin(t),math.cos(t)+math.sin(t)])


def linear_sim(index, t, x, u, dt=1e-4):

    j, t_, x_, u_ = 0, [0], [x], [u]

    while j*dt < t:

        t_.append((j+1)*dt)
        u_.append(get_u(j*dt, index))

        x1 = x_[-1][0] + dt*(x_[-1][2])
        x2 = x_[-1][1] + dt*(x_[-1][3])
        v1 = x_[-1][2] + dt*(-K*x_[-1][2] + math.sqrt(2)*u_[-1][0])
        v2 = x_[-1][3] + dt*(-K*x_[-1][3] + math.sqrt(2)*u_[-1][0])
        theta = x_[-1][4] + dt*(0 + u_[-1][1])

        x_.append(np.array([x1,x2,v1,v2,theta]))

        j += 1
    return np.array(t_),x_


def f(t,x,u):
    k = K
    x1,x2,v1,v2,theta = x
    u1,u2 = u
    return np.array([v1,v2,u1*math.cos(theta)-k*v1,u1*math.sin(theta)-k*v2,u2])

def h(t,x,u):
    x1,x2,v1,v2,theta = x
    u1,u2 = u
    return np.array([x1,x2,theta])

def nonlinear_sim(index,t, x, u, dt=1e-4):
    j, t_, x_,u_ = 0, [0], [x],[u]
    while j*dt < t:
        t_.append((j+1)*dt)
        u_.append(get_u(j*dt, index))
        x1,x2,v1,v2,theta = x_[-1]+dt*f(j*dt,x_[-1],u_[-1])
        x_.append(np.array([x1,x2,v1,v2,theta]))
        j += 1
    return np.array(t_),np.array(x_)




index = 1
t_,x1_ = linear_sim(index,total_t,Xo_1,uo,dt=4e-2)
plt.figure()
lineObjects = plt.plot(t_,x1_,'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.legend(iter(lineObjects),('x1','x2','v1','v2','theta'))
plt.title('Linear State-Time plot 1')


index = 2
t_,x2_ = linear_sim(index,total_t,Xo_2,uo,dt=4e-2)
plt.figure()
lineObjects = plt.plot(t_,x2_,'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.legend(iter(lineObjects),('x1','x2','v1','v2','theta'))
plt.title('Linear State-Time plot 2')


plt.figure()
lineObjects = plt.plot(t_,np.add(x1_,x2_),'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.legend(iter(lineObjects),('x1','x2','v1','v2','theta'))
plt.title('Linear State-Time plot : traj(x1)+traj(x2)')


index = 3
t_,x3_ = linear_sim(index,total_t,np.add(Xo_1,Xo_2),uo,dt=4e-2)
plt.figure()
lineObjects = plt.plot(t_,x3_,'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.legend(iter(lineObjects),('x1','x2','v1','v2','theta'))
plt.title('Linear State-Time plot : traj(x1+x2)')

#------------------------------------------------
index = 1
t_,x1_ = nonlinear_sim(index,total_t,np.add(Xo_1,Xeq),uo,dt=4e-2)
plt.figure()
lineObjects = plt.plot(t_,x1_,'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.legend(iter(lineObjects),('x1','x2','v1','v2','theta'))
plt.title('Nonlinear State-Time plot 1')


index = 2
t_,x2_ = nonlinear_sim(index,total_t,np.add(Xeq,Xo_2),uo,dt=4e-2)
plt.figure()
lineObjects = plt.plot(t_,x2_,'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.legend(iter(lineObjects),('x1','x2','v1','v2','theta'))
plt.title('Nonlinear State-Time plot 2')


plt.figure()
lineObjects = plt.plot(t_,np.add(x1_,x2_),'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.legend(iter(lineObjects),('x1','x2','v1','v2','theta'))
plt.title('Nonlinear State-Time plot : traj(x1)+traj(x2)')


index = 3
t_,x3_ = nonlinear_sim(index,total_t,np.add(np.add(Xo_1,Xo_2),np.add(Xeq,Xeq)),uo,dt=4e-2)
plt.figure()
lineObjects = plt.plot(t_,x3_,'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.legend(iter(lineObjects),('x1','x2','v1','v2','theta'))
plt.title('Nonlinear State-Time plot : traj(x1+x2)')

plt.show()
