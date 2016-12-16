from numpy.linalg import matrix_rank
from numpy.linalg import pinv
from scipy.linalg import expm as sp_expm
import numpy as np
import pylab as plt
import control
import math
import random
import copy

Xo = [[0.],[45.*math.pi/180],[0.]]
Xf = [[0.],[0.],[0.]]
total_time = 5
ts = 5e-5
n_ts = 20000
n_inputs = 2
n_outputs = 3
n_states = 3
K = 1
wait_a_while = 10000
PHYSICAL_DISTURB = 0
OBSERVE_ERROR = 0.001

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

def sim_hat(f, t, x, x_hat, u, K_matrix, L_matrix, dt=ts):
    j_sim, t_sim_, x_sim_, x_hat_sim_, Dy_sim_ = 0, [0], [x], [x_hat], [x_hat-x]
    while j_sim*dt < t:
	# update time
        t_sim_.append((j_sim+1)*dt)
	# compute input
        u1,u2 = -np.dot(K_matrix,x_hat_sim_[-1])
	u = np.array([u1,u2])
        x1,x2,x3 = x_sim_[-1]+dt*f(j_sim*dt,x_sim_[-1],u)
	if x2 > math.pi:
		x2 = x2-2*math.pi
	elif x2 < -math.pi:
		x2 = x2+2*math.pi
        x_sim_.append(np.array([x1,x2,x3]))
	Dy1,Dy2,Dy3 = Dy_sim_[-1]
	x_hat_next = np.dot(Ad,x_hat_sim_[-1])+np.dot(Bd,u)+np.dot(L_matrix,[Dy1,Dy2,Dy3])
	x_hat_sim_.append([x_hat_next[0,0],x_hat_next[1,0],x_hat_next[2,0]])
	Dy_sim_.append(Dy_sim_[-1]+np.dot(A-np.dot(L_matrix,C),Dy_sim_[-1]-np.dot(C,[x1,x2,x3])))
        j_sim += 1
    return np.array(t_sim_),np.array(x_sim_),np.array(x_hat_sim_),np.array(Dy_sim_)


def sim(f, t, x, u, K_matrix, dt=ts):
    j_sim, t_sim_, x_sim_, u_sim_ = 0, [0], [x], [u]
    while j_sim*dt < t:
        t_sim_.append((j_sim+1)*dt)
        u1,u2 = -np.dot(K_matrix,x_sim_[-1])
	u_sim_.append(np.array([u1,u2]))
        x1,x2,x3 = x_sim_[-1]+dt*f(j_sim*dt,x_sim_[-1],u_sim_[-1])
        x_sim_.append(np.array([x1,x2,x3]))
        j_sim += 1
	if j_sim % wait_a_while == 1:
	    h(t,[x1,x2,x3],u_sim_[-1])
    return np.array(t_sim_),np.array(x_sim_),np.array(u_sim_)


def physical_disturbance(x_next_):
	tmp = x_next_
	for i in range(0,n_states):
		tmp[i] = x_next_[i] + random.random()*2*PHYSICAL_DISTURB - PHYSICAL_DISTURB
	return tmp

def observation_error(x_hat_next_):
	tmp = x_hat_next_
	for i in range(0,n_states):
		tmp[i] = x_hat_next_[i] + random.random()*2*OBSERVE_ERROR - OBSERVE_ERROR
	return tmp

# try "iterations" steps to reach the goal:
iterations = int(total_time*n_ts)


A = [[-5.65368219902984, 28.7017492282816, -0.000166867558430104], 
     [                0,                0,                     1], 
     [-8.34425988534470, 132.510363075114, -0.000770248419814788]]

B = [[113.072321034838, 0], 
     [               0, 0], 
     [166.883433779215, 0]]

C = [[1.,0.,0.],
     [0.,1.,0.],
     [0.,0.,1.]]

D = [[  0., 0.],
     [  0., 0.],
     [  0., 0.]]

zero_padder = np.array([[0,0,0,0,0],[0,0,0,0,0]])  
AB = np.concatenate((A,B), axis = 1)
mat = np.concatenate((AB,zero_padder), axis = 0)
exmat = sp_expm(ts*mat)

Ad = exmat[0:n_states,0:n_states]
Bd = exmat[0:n_states,n_states:n_states+n_inputs]

print "\n -------------Closed-Loop Control------------- "

eigval,eigvec = np.linalg.eig(A)
print "\neigen values of A = " + str(eigval)

# -12.87243308  -3.64492168  10.86290231
K_matrix = control.place(A,B,[-12.87243308,-3.64492168,-10.86290231])
print "\nK matrix = \n" + str(K_matrix)

A_BK = A - np.dot(B,K_matrix)
print "\nA - BK = \n" + str(A_BK)

eigval_newC,eigvec_newC = np.linalg.eig(A_BK)
print "\neigen values of (A-BK) = " + str(eigval_newC)

print "\n now lets do some simulations:"
t_ = np.zeros((1,iterations))
x_ = np.zeros((n_states,iterations+1))
z_ = np.zeros((1,iterations+1))
u_ = np.zeros((n_inputs,iterations+1))

x0 = np.array([[Xo[0][0]],[Xo[1][0]],[Xo[2][0]]])
x_[:,0] = x0[:,0]

for i in range(0,iterations):
	t_[:,i] = ts*i
	u_[:,i+1] = -np.dot(K_matrix,[x_[0,i],x_[1,i],x_[2,i]])
	x_next = np.dot(Ad,x_[:,i])+np.dot(Bd,u_[:,i+1])
	if x_next[1] > math.pi:
		x_next[1] = x_next[1]-2*math.pi
	elif x_next[1] < -math.pi:
		x_next[1] = x_next[1]+2*math.pi
	x_[:,i+1] = x_next
t_ = np.hstack((np.squeeze(t_),np.array([iterations*ts])))

# linear simulation
plt.figure()
plt.plot(t_, x_[0,:])
plt.plot(t_, x_[1,:])
plt.plot(t_, x_[2,:])
plt.plot(t_, u_[0,:])
plt.plot(t_, u_[1,:])
plt.legend(['x1','x2','x3','u1','u2'],loc=4)
plt.title('Linear simulation')
plt.plot(t_,z_[0,:],'k--')

# nonlinear simulation
plt.figure()
x0 = np.array([Xo[0][0],Xo[1][0],Xo[2][0]])
u0 = np.array([0,0])
t_sim_,x_sim_,u_sim_ = sim(f,total_time,x0,u0,K_matrix,dt=ts)

plt.plot(t_sim_,x_sim_,'-')
plt.plot(t_sim_,u_sim_,'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.legend(['x1','x2','x3','u1','u2'],loc=4)
plt.plot(t_sim_,z_[0,:],'k--')
plt.title('Nonlinear simulation')

print "\n now lets compare linear and nonlinear results:"
x_difference_ = x_sim_-np.transpose(x_)
u_difference_ = u_sim_-np.transpose(u_)

# differnce plot
t_log_ = [t_sim_[0]]
for i in range(1,iterations+1):
	t_log_.append(math.log(t_sim_[i]))

plt.figure()
plt.plot(t_log_,x_difference_,'-')
plt.plot(t_log_,u_difference_,'-')
plt.legend(['x1','x2','x3','u1','u2'],loc = 4)
plt.xlabel('time (log (sec))')
plt.ylabel('difference magnitude')
plt.title('Difference between Linear and Nonlinear Simulation (log scaled time)')
plt.plot(t_,z_[0,:],'k--')

plt.figure()
plt.plot(t_sim_,x_difference_,'-')
plt.plot(t_sim_,u_difference_,'-')
plt.legend(['x1','x2','x3','u1','u2'],loc = 4)
plt.xlabel('time (sec)')
plt.ylabel('difference magnitude')
plt.title('Difference between Linear and Nonlinear Simulation')
plt.plot(t_,z_[0,:],'k--')


print "\n -------------Closed-Loop Observation------------- "

eigval,eigvec = np.linalg.eig(A)
print "\neigen values of A = " + str(eigval)

# [-12.87243308  -3.64492168  10.86290231]
L_T = control.place(np.transpose(A),np.transpose(C),[-12.87243308,-3.64492168 ,-10.86290231])
L_matrix = np.transpose(L_T)

L_matrix = [[-3.8537, 28.7017,-0.0002], # eigen values: [-1.87243308,-1.64492168,-1.86290231]
            [      0,     1.6,    1.0],
            [-8.3443,132.5104, 1.8592]]
'''
L_matrix = [[7.2187,  28.7017, -0.0002], # eigen values:[-12.87243308, -3.64492168, -10.86290231]
            [      0,  3.6449,     1.0],
            [-8.3443,132.5104, 10.8621]]
'''
print "\nL matrix = \n" + str(L_matrix)

A_LC = A - np.dot(L_matrix,C)
print "\nA - LC = \n" + str(A_LC)

eigval_newO,eigvec_newO = np.linalg.eig(A_LC)
print "\neigen values of (A_LC) = " + str(eigval_newO)

print "\n now lets do some simulations:"
# plot fixed error
t_ = np.zeros((1,iterations))
error_ = np.zeros((n_states,iterations+1))

x_real = np.array([[Xo[0][0]],[Xo[1][0]],[Xo[2][0]]])
x_hat = np.array([[Xo[0][0]+20],[Xo[1][0]+1],[Xo[2][0]+6]])
error_[:,0] = x_real[:,0] - x_hat[:,0]

print "xo_real = " + str(x_real)
print "xo_real = " + str(x_hat)
print "xo_error = " + str(error_[:,0])

for i in range(0,iterations):
	t_[:,i] = ts*i
	Derror = np.dot(A_LC,error_[:,i])
	if i < 100 and i%10 == 1:
		print "\nerror = " + str(error_[:,i])
		print "Derror = " + str(Derror)
		print "next_error  = " + str(error_[:,i]+Derror)
	error_[:,i+1] = error_[:,i]+Derror
t_ = np.hstack((np.squeeze(t_),np.array([iterations*ts])))

plt.figure()
plt.plot(t_, error_[0,:])
plt.plot(t_, error_[1,:])
plt.plot(t_, error_[2,:])
plt.legend(['error x1','error x2','error x3'],loc=4)
plt.xlim([0,0.002])
plt.xlabel('time (sec)')
plt.ylabel('error')
plt.title('The Error Plot (catching up on fixed error)')
plt.plot(t_,z_[0,:],'k--')

# plot free motion error
tmp = iterations
iterations = 3700
t_ = np.zeros((1,iterations))
x_real_ = np.zeros((n_states,iterations+1))
x_hat_ = np.zeros((n_states,iterations+1))
error_ = np.zeros((n_states,iterations+1))
z_small_ = np.zeros((1,iterations+1))

x_real_[:,0] = x_real[:,0]
x_hat_[:,0] = x_hat[:,0]
error_[:,0] = x_hat[:,0]-x_real[:,0]

x_hat_next = x_hat_[:,0]
x_next = x_real_[:,0]

for i in range(0,iterations):
	t_[:,i] = ts*i
	x_next = np.dot(Ad,x_real_[:,i])
	x_hat_next = np.dot(Ad,x_hat_[:,i])-0.001*(x_hat_next-x_next)
	#x_hat_next = np.dot(Ad,x_hat_[:,i])-np.dot(L_matrix,np.dot(C,error_[:,i]))

	if x_next[1] > math.pi:
		x_next[1] = x_next[1]-2*math.pi
	elif x_next[1] < -math.pi:
		x_next[1] = x_next[1]+2*math.pi

	if x_hat_next[1] > math.pi:
		x_hat_next[1] = x_hat_next[1]-2*math.pi
	elif x_hat_next[1] < -math.pi:
		x_hat_next[1] = x_hat_next[1]+2*math.pi
	
	Derror = np.dot(A-np.dot(L_matrix,C),error_[:,i])
	error_[:,i+1] = x_hat_next-x_next
	#error_[:,i+1] = error_[:,i]+Derror
	x_hat_[:,i+1] = x_hat_next
	x_real_[:,i+1] = x_next

for i in range(0,iterations+1):
	if x_real_[1,i] < 0:
		x_real_[1,i] = x_real_[1,i]+2*math.pi
	if x_hat_[1,i] < 0:
		x_hat_[1,i] = x_hat_[1,i]+2*math.pi

t_ = np.hstack((np.squeeze(t_),np.array([iterations*ts])))

f_, axarr = plt.subplots(3, sharex = True)
axarr[0].plot(t_, x_real_[0,:])
axarr[0].plot(t_, x_hat_[0,:])
axarr[0].plot(t_, error_[0,:])
axarr[0].legend(['x1','x_hat1','error x1'],loc=1)
axarr[0].set_title('The Error Plot (catching up on initial error)')
axarr[0].plot(t_, z_small_[0,:],'--k')

axarr[1].plot(t_, x_real_[1,:])
axarr[1].plot(t_, x_hat_[1,:])
axarr[1].plot(t_, error_[1,:])
axarr[1].legend(['x2','x_hat2','error x2'],loc=4)
axarr[1].set_title('The Error Plot (catching up on initial error)')
axarr[1].plot(t_, z_small_[0,:],'--k')

axarr[2].plot(t_, x_real_[2,:])
axarr[2].plot(t_, x_hat_[2,:])
axarr[2].plot(t_, error_[2,:])
axarr[2].legend(['x3','x_hat3','error x3'],loc=4)
axarr[2].set_title('The Error Plot (catching up on initial error)')
axarr[2].plot(t_, z_small_[0,:],'--k')

iterations = tmp

print "\n -------------Closed-Loop Everything----------------- "
zero_padder = np.zeros((3,3))
left = np.concatenate((A_LC,np.dot(L_matrix,C)),axis = 0)
right = np.concatenate((zero_padder,A_BK),axis = 0)

Big_Matrix = np.concatenate((left,right),axis = 1)
print "Big_Matrix = \n" + str(Big_Matrix)

eigval_new,eigvec_new = np.linalg.eig(Big_Matrix)
print "\neigen values of Big_Matrix = \n" + str(eigval_new)


print "\n now lets do some simulations:"  
t_ = np.zeros((1,iterations))
x_ = np.zeros((n_states,iterations+1))
x_hat_ = np.zeros((n_states,iterations+1))
Dy_ = np.zeros((n_states,iterations+1))
z_ = np.zeros((1,iterations+1))
u_ = np.zeros((n_inputs,iterations+1))

x0 = np.array([[Xo[0][0]],[Xo[1][0]],[Xo[2][0]]])
x0_hat = np.array([[Xo[0][0]-26.0],[Xo[1][0]+1.0],[Xo[2][0]+23.0]])
x_[:,0] = x0[:,0]
x_hat_[:,0] = x0_hat[:,0]
Dy_[:,0] = x_[:,0] - x_hat_[:,0]

for i in range(0,iterations):
	t_[:,i] = ts*i
	u_[:,i+1] = -np.dot(K_matrix,[x_hat_[0,i],x_hat_[1,i],x_hat_[2,i]])
	x_next = np.dot(Ad,x_[:,i])+np.dot(Bd,u_[:,i+1])
	x_hat_next = np.dot(Ad,x_hat_[:,i])+np.dot(Bd,u_[:,i+1])+np.dot(L_matrix,Dy_[:,i])
	x_[:,i+1] = x_next
	x_hat_[:,i+1] = x_hat_next
	Dy_[:,i+1] = Dy_[:,i]+np.dot(A-np.dot(L_matrix,C),error_[:,i])

t_ = np.hstack((np.squeeze(t_),np.array([iterations*ts])))

# linear simulation
ff, axarrr = plt.subplots(2, sharex = True)
axarrr[0].plot(t_, x_[0,:])
axarrr[0].plot(t_, x_[1,:])
axarrr[0].plot(t_, x_[2,:])
axarrr[0].plot(t_, x_hat_[0,:])
axarrr[0].plot(t_, x_hat_[1,:])
axarrr[0].plot(t_, x_hat_[2,:])
axarrr[0].legend(['x1','x2','x3','x_hat1','x_hat2','x_hat3'],loc=4)
axarrr[0].set_title('Controller/Observer State Plot (Linear)')
axarrr[0].plot(t_,z_[0,:],'k--')

axarrr[1].plot(t_, Dy_[0,:])
axarrr[1].plot(t_, Dy_[1,:])
axarrr[1].plot(t_, Dy_[2,:])
axarrr[1].legend(['Dy1','Dy2','Dy3'],loc=4)
axarrr[1].set_title('Controller/Observer Error Plot (Linear)')
axarrr[1].plot(t_,z_[0,:],'k--')

# nonlinear simulation
x0 = np.array([[Xo[0][0]],[Xo[1][0]],[Xo[2][0]]])
x0_hat = np.array([[Xo[0][0]+0.1],[Xo[1][0]+0.1],[Xo[2][0]-0.1]])
u0 = np.array([0,0])
t_sim_,x_sim_,x_hat_sim_,Dy_sim_ = sim_hat(f,total_time/5,x0,x0_hat,u0,K_matrix,L_matrix,dt=ts)

fff, axarrrr = plt.subplots(2, sharex = True)
np.reshape(x_sim_,(x_sim_.shape[0], 3))
axarrrr[0].plot(t_sim_,x_sim_[:,:,0],'-')
axarrrr[0].plot(t_sim_,x_hat_sim_,'-')
axarrrr[0].legend(['x1','x2','x3','x_hat1','x_hat2','x_hat3'],loc=4)
axarrrr[0].set_title('Controller/Observer State Plot (Nonlinear)')
#axarrrr[0].plot(t_sim_,z_[0,:],'k--')

axarrrr[1].plot(t_sim_,Dy_sim_[:,:,0],'-')
axarrrr[1].legend(['Dy1','Dy2','Dy3'],loc=4)
axarrrr[1].set_title('Controller/Observer Error Plot (Nonlinear)')
#axarrrr[1].plot(t_sim_,z_[0,:],'k--')




plt.show()
