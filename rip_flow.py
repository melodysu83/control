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

def get_u(t):
    u1 = 0.2*math.sin(0.5*t)
    u2 = 0
    return np.array([u1,u2])

def sim(f, t, x, u, dt=1e-4):
    j, t_, x_, u_ = 0, [0], [x], [u]
    while j*dt < t:
        t_.append((j+1)*dt)
        u_.append(get_u(j*dt))
        x1,x2,x3 = x_[-1]+dt*f(j*dt,x_[-1],u_[-1])
	'''
        x2 = x2 % (2*math.pi)
        if x2 > math.pi:
            x2 -= 2*math.pi
	'''
        x_.append(np.array([x1,x2,x3]))
        j += 1
    return np.array(t_),np.array(x_)


# simulation time
total_t = 25

# start state
xs = np.array([0,math.pi/4.0,0])
us = np.array([0,0])

#do the simulation
t_,x_ = sim(f,total_t,xs,us,dt=1e-4)

# plot the state-time plot
plt.figure()
lineObjects = plt.plot(t_,x_,'-')
plt.xlabel('time (sec)')
plt.ylabel('state (init =[0,math.pi/4.0,0])')
plt.legend(iter(lineObjects),('x1 = dAlpha','x2 = Beta','x3 = dBeta'))
plt.title('State-Time plot')

# plot state with different initial condition
xs = np.array([0,0,0])
us = np.array([0,0])
total_t = 15
plt.figure()
for i in range(5):
    xinit = xs+np.array([0,10*i*math.pi/180.0,0])
    t_,x_ = sim(f,total_t,xinit,us,dt=1e-3)
    plt.plot(t_,x_[:,0],'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.title('x1-time plot (with different initial condition)')

plt.figure()
for i in range(5):
    xinit = xs+np.array([0,10*i*math.pi/180.0,0])
    t_,x_ = sim(f,total_t,xinit,us,dt=1e-3)
    plt.plot(t_,x_[:,1],'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.title('x2-time plot (with different initial condition)')

plt.figure()
for i in range(5):
    xinit = xs+np.array([0,10*i*math.pi/180.0,0])
    t_,x_ = sim(f,total_t,xinit,us,dt=1e-3)
    plt.plot(t_,x_[:,2],'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.title('x3-time plot (with different initial condition)')


xs = np.array([0,0,0])
us = np.array([0,0])
total_t = 3
plt.figure()
for i in range(5):
    xinit = xs+np.array([0,10*i*math.pi/180.0,0])
    t_,x_ = sim(f,total_t,xinit,us,dt=1e-3)
    plt.plot(t_,x_[:,0],'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.title('x1-time plot (with different initial condition)')

plt.figure()
for i in range(5):
    xinit = xs+np.array([0,10*i*math.pi/180.0,0])
    t_,x_ = sim(f,total_t,xinit,us,dt=1e-3)
    plt.plot(t_,x_[:,1],'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.title('x2-time plot (with different initial condition)')

plt.figure()
for i in range(5):
    xinit = xs+np.array([0,10*i*math.pi/180.0,0])
    t_,x_ = sim(f,total_t,xinit,us,dt=1e-3)
    plt.plot(t_,x_[:,2],'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.title('x3-time plot (with different initial condition)')


total_t = 15
xs = [0,0,0]
us = [0,0]
print 'f(x0) = ', f(0.,xs,us)

# part 1:
t_,x_= sim(f,total_t,xs,us,dt=1e-3)
for i in range(1,6):
    xinit = xs+np.array([0,10*i*math.pi/180.0,0])
    t_,x_ = sim(f,total_t,xinit,us,dt=1e-3)
    print "a=" + str(i) + " : x(t)=" , x_[total_t/(1e-3)-1,:] 

# Here is the result for part 1 #
xta1 = np.array([3.51043591,3.12480531,0.01556264])
xta2= np.array( [3.51062021,3.12480658,0.01460275])
xta3= np.array( [3.51067900,3.12475983,0.01410690])
xta4= np.array( [3.51063927,3.12469731,0.01406350])
xta5= np.array( [3.51052790,3.12464651,0.01444228])

# part 2:
a_ = np.array([1,2,3,4])
for i in range(3):
    plt.figure()
    lineObjects = plt.plot(a_,np.array([(xta2[i]-xta1[i])/1.0,(xta3[i]-xta2[i])/2.0,(xta4[i]-xta3[i])/3.0,(xta5[i]-xta4[i])/4.0]),'.-')
    plt.title('x'+str(i+1))

plt.show()
 
