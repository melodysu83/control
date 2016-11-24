import numpy as np
import pylab as plt
import math
import copy

K = 1

A = [[28.7009,-5.635,-0.0002],[0,0,1],[132.5065,-8.3440,-0.0008]]
B = [[113.0690,0],[0,0],[166.8785,1]]

total_t = 0.2 
total_t_nonlinear = 3 

#               dAlpha,         Beta,   dBeta
Xeq  = np.array([    0,            0,      0]) 
Xo_1 = np.array([    0, math.pi/18.0,      0]) 
Xo_2 = np.array([0.001,            0,      0]) 

ueq = np.array([0,0])
uo  = np.array([0,0])

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
        u1 = math.sin(0.001*t)
        u2 = 0
        return np.array([u1,u2])
    elif index == 2:
        u1 = math.cos(0.001*t)
        u2 = 0
        return np.array([u1,u2])
    else:
        return np.array([math.cos(0.001*t)+math.sin(0.001*t),0])


def linear_sim(index, t, x, u, dt=1e-4):

    j, t_, x_, u_ = 0, [0], [x], [u]

    while j*dt < t:

        t_.append((j+1)*dt)
        u_.append(get_u(j*dt, index))
        x1 = x_[-1][0] + dt*(A[0][0]*x_[-1][0]+A[0][1]*x_[-1][1]+A[0][2]*x_[-1][2]+B[0][0]*u_[-1][0]+B[0][1]*u_[-1][1])
        x2 = x_[-1][1] + dt*(A[1][0]*x_[-1][0]+A[1][1]*x_[-1][1]+A[1][2]*x_[-1][2]+B[1][0]*u_[-1][0]+B[1][1]*u_[-1][1])
        x3 = x_[-1][2] + dt*(A[2][0]*x_[-1][0]+A[2][1]*x_[-1][1]+A[2][2]*x_[-1][2]+B[2][0]*u_[-1][0]+B[2][1]*u_[-1][1])
        x_.append(np.array([x1,x2,x3]))
        j += 1
    return np.array(t_),np.array(x_)


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


def nonlinear_sim(index,t, x, u, dt=1e-4):
    j, t_, x_, u_ = 0, [0], [x], [u]
    while j*dt < t:
        t_.append((j+1)*dt)
        u_.append(get_u(j*dt,index))
        x1,x2,x3 = x_[-1]+dt*f(j*dt,x_[-1],u_[-1])
        x_.append(np.array([x1,x2,x3]))
        j += 1
    return np.array(t_),np.array(x_)




index = 1
t_,x1_ = linear_sim(index,total_t,Xo_1,uo,dt=1e-4)
plt.figure()
lineObjects = plt.plot(t_,x1_,'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.legend(iter(lineObjects),('x1 = dAlpha','x2 = Beta','x3 = dBeta'),loc='lower left')
plt.title('Linear State-Time plot 1')


index = 2
t_,x2_ = linear_sim(index,total_t,Xo_2,uo,dt=1e-4)
plt.figure()
lineObjects = plt.plot(t_,x2_,'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.legend(iter(lineObjects),('x1 = dAlpha','x2 = Beta','x3 = dBeta'),loc='lower left')
plt.title('Linear State-Time plot 2')


plt.figure()
lineObjects = plt.plot(t_,np.add(x1_,x2_),'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.legend(iter(lineObjects),('x1 = dAlpha','x2 = Beta','x3 = dBeta'),loc='lower left')
plt.title('Linear State-Time plot : traj(x1)+traj(x2)')


index = 3
t_,x3_ = linear_sim(index,total_t,np.add(Xo_1,Xo_2),uo,dt=1e-4)
plt.figure()
lineObjects = plt.plot(t_,x3_,'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.legend(iter(lineObjects),('x1 = dAlpha','x2 = Beta','x3 = dBeta'),loc='lower left')
plt.title('Linear State-Time plot : traj(x1+x2)')

#------------------------------------------------
index = 1
t_,x1_ = nonlinear_sim(index,total_t_nonlinear,np.add(Xo_1,Xeq),uo,dt=1e-4)
plt.figure()
lineObjects = plt.plot(t_,x1_,'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.legend(iter(lineObjects),('x1 = dAlpha','x2 = Beta','x3 = dBeta'))
plt.title('Nonlinear State-Time plot 1')


index = 2
t_,x2_ = nonlinear_sim(index,total_t_nonlinear,np.add(Xeq,Xo_2),uo,dt=1e-4)
plt.figure()
lineObjects = plt.plot(t_,x2_,'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.legend(iter(lineObjects),('x1 = dAlpha','x2 = Beta','x3 = dBeta'))
plt.title('Nonlinear State-Time plot 2')


plt.figure()
lineObjects = plt.plot(t_,np.add(x1_,x2_),'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.legend(iter(lineObjects),('x1 = dAlpha','x2 = Beta','x3 = dBeta'))
plt.title('Nonlinear State-Time plot : traj(x1)+traj(x2)')


index = 3
t_,x3_ = nonlinear_sim(index,total_t_nonlinear,np.add(np.add(Xo_1,Xo_2),np.add(Xeq,Xeq)),uo,dt=1e-4)
plt.figure()
lineObjects = plt.plot(t_,x3_,'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.legend(iter(lineObjects),('x1 = dAlpha','x2 = Beta','x3 = dBeta'))
plt.title('Nonlinear State-Time plot : traj(x1+x2)')

plt.show()
