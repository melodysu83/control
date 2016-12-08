from numpy.linalg import matrix_rank
from numpy.linalg import pinv
from scipy.linalg import expm as sp_expm
import numpy as np
import pylab as plt
import math
import svd

total_time = 0.1
ts = 5e-5
n_ts = 20000
n_inputs = 2
n_states = 3
time_for_reaching_subgoal = 100
K = 1

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
    x1,x2,x3 = x
    u1,u2 = u
    print "h(t,x,u) = " + str(np.array([x1,x2,x3]))
    return np.array([x1,x2,x3])

def sim(f, t, x, u, dt=ts):
    j_sim, t_sim_, x_sim_, u_sim_ = 0, [0], [x], u
    while j_sim*dt < t:
        t_sim_.append((j_sim+1)*dt)
        u_sim_now = np.array([u[n_inputs*j_sim,0],u[n_inputs*j_sim+1,0]])
        x1,x2,x3 = x_sim_[-1]+dt*f(j_sim*dt,x_sim_[-1],u_sim_now)
        x_sim_.append(np.array([x1,x2,x3]))
        j_sim += 1
	if j_sim % time_for_reaching_subgoal == 0:
	    h(t,[x1,x2,x3],u_sim_now)
    return np.array(t_sim_),np.array(x_sim_)



 
A = [[ 28.7009,  -5.6535,  -0.0002],
     [      0.,       0.,       1.],
     [132.5065,  -8.3440,  -0.0008]]

B = [[       113.0690, 0.],
     [             0., 0.],
     [       166.8785, 1.]]

C = [[1.,0.,0.],
     [0.,1.,0.],
     [0.,0.,1.]]

D = [[  0., 0.],
     [  0., 0.],
     [  0., 0.]]


# discretizing

# [[n_states+n_inputs],...n_inputs times...,[n_states+n_inputs]]
zero_padder = np.array([[0,0,0,0,0],[0,0,0,0,0]])  
AB = np.concatenate((A,B), axis = 1)
mat = np.concatenate((AB,zero_padder), axis = 0)
exmat = sp_expm(ts*mat)

Ad = exmat[0:n_states,0:n_states]
Bd = exmat[0:n_states,n_states:n_states+n_inputs]


print "(1) Controllability check:"
# try "iterations" steps to reach the goal:
iterations = int(total_time*n_ts)

Control_matrix = np.zeros((n_states,iterations*n_inputs))

for i in range(0,iterations):
	Control_matrix[:,i*n_inputs:(i+1)*n_inputs] = np.dot((np.linalg.matrix_power(Ad,iterations-i)),Bd)

Control_rank = matrix_rank(Control_matrix)
print "rank = " + str(Control_rank)  # rank = 3 --> fully controllable.


print "(2) find steering input sequence"
Control_matrix_svd = np.linalg.svd(Control_matrix)

zero_padder = np.zeros((n_states,iterations*n_inputs-n_states))

V = Control_matrix_svd[2]
S = np.concatenate((np.diag(Control_matrix_svd[1]), zero_padder),axis = 1)
U = Control_matrix_svd[0]

# psudo inverse
s = np.linalg.pinv(S)

Xo = [[0.],[5.*math.pi/180],[0.]]
Xf = [[0.],[0.],[0.]]

Z = Xf - np.dot(np.linalg.matrix_power(Ad,iterations),Xo)
print "Z = " + str(Z)


# this is our control input
u = np.dot(np.dot(np.dot(np.transpose(V),s),np.transpose(U)),Z)



print "(3) See the simulation result:"
t_ = np.zeros((1,iterations))
x_ = np.zeros((n_states,iterations+1))
z_ = np.zeros((1,iterations+1))

x0 = np.array([[Xo[0][0]],[Xo[1][0]],[Xo[2][0]]])
x_[:,0] = x0[:,0]

for i in range(0,iterations):
	t_[:,i] = ts*i
	x_[:,i+1] = np.dot(Ad,x_[:,i])+np.squeeze(np.dot(Bd,u[n_inputs*i:n_inputs*(i+1),:]))
t_ = np.hstack((np.squeeze(t_),np.array([iterations*ts])))

plt.figure()
plt.plot(t_, x_[0,:])
plt.plot(t_, x_[1,:])
plt.plot(t_, x_[2,:])
plt.legend(['x1','x2','x3'],loc=2)
plt.title('Linear simulation')
plt.plot(t_,z_[0,:],'k--')




plt.figure()
x0 = np.array([Xo[0][0],Xo[1][0],Xo[2][0]])
t_,x_ = sim(f,total_time,x0,u,dt=ts)

lineObjects = plt.plot(t_,x_,'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.legend(iter(lineObjects),('x1','x2','x3'),loc=2)
plt.plot(t_,z_[0,:],'k--')
plt.title('Nonlinear simulation')
plt.show()
