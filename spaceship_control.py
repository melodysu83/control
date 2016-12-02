from numpy.linalg import matrix_rank
from numpy.linalg import pinv
import numpy as np
import pylab as plt
import math
import svd

def matmul(mat,vec,lenth):
	matvec = [[vec[0][0],0],[vec[1][0],0],[vec[2][0],0],[vec[3][0],0],[vec[4][0],0]]
	result = svd.matrixmultiply(mat,matvec)
	if lenth == 5:
		return [[result[0][0]],[result[1][0]],[result[2][0]],[result[3][0]],[result[4][0]]]
	elif lenth == 10:
		return [[result[0][0]],[result[1][0]],[result[2][0]],[result[3][0]],[result[4][0]],
			[result[5][0]],[result[6][0]],[result[7][0]],[result[8][0]],[result[9][0]]]

def f(t,x,u):
    x1,x2,v1,v2,theta = x
    u1,u2 = u
    return np.array([v1,v2,u1*math.cos(theta)-K*v1,u1*math.sin(theta)-K*v2,u2])

def h(t,x,u):
    x1,x2,v1,v2,theta = x
    u1,u2 = u
    print "h(t,x,u) = " + str(np.array([x1,x2,theta]))
    return np.array([x1,x2,theta])

def get_u(j):
    u1 = U[2*(4-j)][0]
    u2 = U[2*(4-j)+1][0]
    return np.array([u1,u2])

def sim(f, t, x, u, dt=1e-4):
    j, t_, x_, u_ = 0, [0], [x], [u]
    while j*dt < t:
        t_.append((j+1)*dt)
	u_now = get_u(int(j/time_for_reaching_subgoal))
        u_.append(u_now)
        x1,x2,v1,v2,theta = x_[-1]+dt*f(j*dt,x_[-1],u_[-1])
        theta = theta % (2*math.pi)
        if theta > math.pi:
            theta -= 2*math.pi
        x_.append(np.array([x1,x2,v1,v2,theta]))
        j += 1
	if j % time_for_reaching_subgoal == 0:
	    h(t,[x1,x2,v1,v2,theta],[u_now[0],u_now[1]])
    return np.array(t_),np.array(x_)

K = 1
zero = [0,0,0,0,0]
time_for_reaching_subgoal = 40

A = [[0.,0.,1., 0.,              0.],
     [0.,0.,0., 1.,              0.],
     [0.,0.,-K, 0.,-1./math.sqrt(2)],
     [0.,0.,0., -K, 1./math.sqrt(2)],
     [0.,0.,0., 0.,              0.]]

B = [[             0., 0.],
     [             0., 0.],
     [1./math.sqrt(2), 0.],
     [1./math.sqrt(2), 0.],
     [             0., 1.]]

C = [[1.,0.,0.,0.,0.],
     [0.,1.,0.,0.,0.],
     [0.,0.,0.,0.,1.]]

D = [[        0., 0.],
     [        0., 0.],
     [        0., 0.]]

# initial condition:  Xo = [1,1,0,0,pi/4]
#                     Uo = [1,1]
# final condition:    Xf = [3,8,0,0,pi]
#                     Uf = [0,0]

print "(1) Controllability check:"
AB = svd.matrixmultiply(A,B)
A2B = svd.matrixmultiply(A,AB)
A3B = svd.matrixmultiply(A,A2B)
A4B = svd.matrixmultiply(A,A3B)

Control_matrix = np.concatenate((A4B,A3B,A2B,AB,B),axis=1)
print Control_matrix
# This is the controllability matrix I get:
#
# [[-0.70710678 -0.70710678  0.70710678  0.70710678 -0.70710678 -0.70710678   0.70710678  0.          0.          0.        ]
#  [-0.70710678  0.70710678  0.70710678 -0.70710678 -0.70710678  0.70710678   0.70710678  0.          0.          0.        ]
#  [ 0.70710678  0.70710678 -0.70710678 -0.70710678  0.70710678  0.70710678  -0.70710678 -0.70710678  0.70710678  0.        ]
#  [ 0.70710678 -0.70710678 -0.70710678  0.70710678  0.70710678 -0.70710678  -0.70710678  0.70710678  0.70710678  0.        ]
#  [ 0.          0.          0.          0.          0.          0.          0.           0.          0.          1.        ]]



Control_rank = matrix_rank(Control_matrix)
print "rank = " + str(Control_rank)  # rank = 5 --> fully controllable.


print "(2) find steering input sequence"
# we need to transpose first before feeding into svd function
# due to shape constraint
a,b,c = svd.svd(svd.transpose(Control_matrix))

# This is the "a" I get:
#
#[[   -0.48324982438233488,  7.6327832942979512e-16, -1.0178532794652785e-14,     0.12833396757851195,  0.0], 
# [ 6.9388939039072284e-16,     0.55105881968666859,     0.17224259223220284,   1.385436904088877e-14,  0.0], 
# [    0.48324982438233482, -8.6042284408449632e-16,  1.0440476039887174e-14,    -0.12833396757851209,  0.0], 
# [-7.4940054162198066e-16,    -0.55105881968666848,    -0.17224259223220295, -1.3783679059242715e-14,  0.0], 
# [   -0.48324982438233488,  7.4940054162198066e-16, -1.0398842676101839e-14,     0.12833396757851209,  0.0], 
# [ 7.4940054162198066e-16,     0.55105881968666848,     0.17224259223220295,  1.3783679059242715e-14,  0.0], 
# [    0.48324982438233488, -7.4940054162198066e-16,  1.0398842676101839e-14,    -0.12833396757851209,  0.0], 
# [ -5.134781488891349e-16,    -0.29833292097354391,     0.95446187365624646,  7.6645287339083268e-14,  0.0], 
# [   -0.25666793515702424,  4.8572257327350599e-16,  7.7440328390661301e-14,    -0.96649964876466921,  0.0], 
# [                    0.0,                     0.0,                    -0.0,                     0.0, -1.0]]


# This is the "b" I get:
#
#[2.92080962648189, 2.557612414958357, 0.6772139505731477, 0.6847416489820997, 1.0]


# This is the "c" I get:
#
#[[  0.4679650802706302, -0.45705607224241235,  -0.5395366037872631,  -0.5301025218270025, -0.0], 
# [ 0.46796508027063127,  0.45705607224241102,  0.53953660378734702, -0.53010252182691653, -0.0], 
# [-0.53010252182695883,  0.53953660378730584, -0.45705607224237449, -0.46796508027066708, -0.0], 
# [-0.53010252182696027, -0.53953660378730417,  0.45705607224244887, -0.46796508027059447, -0.0], 
# [                 0.0,                  0.0,                 -0.0,                 -0.0, -1.0]]



a = np.concatenate((a,[zero,zero,zero,zero,zero,zero,zero,zero,zero,zero]),axis=1)

b = np.concatenate((np.diag(b),[zero,zero,zero,zero,zero]),axis=0)

v_star = svd.transpose(a)
w = svd.transpose(b)
u = c

print "Now this should be identical to the Control_matrix:"
print svd.matrixmultiply(u,svd.matrixmultiply(w,v_star)) # Control_matrix = u * w * v_star

Xo = [[1.],[1.],[0.],[0.],[math.pi/4.]]
Xf = [[1.15],[1.15],[0.1],[0.1],[1]]

A5Xo = matmul(svd.matrixmultiply(A,svd.matrixmultiply(svd.matrixmultiply(A,A),svd.matrixmultiply(A,A))),Xo,5)

Z = [[Xf[0][0]-A5Xo[0][0]],[Xf[1][0]-A5Xo[1][0]],[Xf[2][0]-A5Xo[2][0]],[Xf[3][0]-A5Xo[3][0]],[Xf[4][0]-A5Xo[4][0]]]

print "Xf = " + str(Xf)
print "A5Xo = " + str(A5Xo)
print "Xf - A5Xo = " + str(Z)

v = pinv(v_star)
u_star = pinv(u)
w_dagger = np.diag(np.reciprocal(np.diag(w)))
w_dagger = np.concatenate((w_dagger,[zero,zero,zero,zero,zero]),axis=0)

U = matmul(svd.matrixmultiply(svd.matrixmultiply(v,w_dagger),u_star),Z,10)

print "U = " + str(U)


print "See the simulation result:"
plt.figure()
total_t = 5*4e-3*time_for_reaching_subgoal
x0 = [Xo[0][0],Xo[1][0],Xo[2][0],Xo[3][0],Xo[4][0]]
u0 = [1,1]
t_,x_ = sim(f,total_t,x0,u0,dt=4e-3)

lineObjects = plt.plot(t_,x_,'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.legend(iter(lineObjects),('x1','x2','v1','v2','theta'),loc=2)

plt.show()

#------------------------------------------------------------

