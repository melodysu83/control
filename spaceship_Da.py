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

# part 1:
t_,x_,x1_,x2_ = sim(f,total_t,x0,u0,dt=4e-2)
for i in range(1,6):
    xinit = x0+np.array([0,0,0,0,10*i*math.pi/180.0])
    t_,x_,x1_,x2_ = sim(f,total_t,xinit,u0,dt=4e-2)
    print "a=" + str(i) + " : x(t)=" , x_[60/(4e-2)-1,:] 

# Here is the result for part 1 #
xta1 = np.array([ 1.82100427,  5.65061528, -0.87326504,  0.46588623,  2.79551686])
xta2= np.array( [ 0.81212008,  5.88098381, -0.94089848,  0.30716749,  2.97004978])
xta3= np.array( [-0.22143997,  5.93266163, -0.97994319,  0.13911562 ,-3.1386026 ])
xta4= np.array( [-1.24827168 , 5.80407852, -0.98921283, -0.0331632,  -2.96406967])
xta5= np.array( [-2.23717529,  5.49914142, -0.96842573, -0.20443439, -2.78953675])

# part 2:
a_ = np.array([1,2,3,4])
for i in range(5):
    plt.figure()
    lineObjects = plt.plot(a_,np.array([(xta2[i]-xta1[i])/1.0,(xta3[i]-xta2[i])/2.0,(xta4[i]-xta3[i])/3.0,(xta5[i]-xta4[i])/4.0]),'.-')
    plt.title('x'+str(i+1))

plt.show()

    

