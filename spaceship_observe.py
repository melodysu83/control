from numpy.linalg import matrix_rank
from numpy.linalg import pinv
import numpy as np
import pylab as plt
import math
import svd

def matmul(mat,vec):
	matvec = [[vec[0][0],0],[vec[1][0],0],[vec[2][0],0],[vec[3][0],0],[vec[4][0],0],
		  [vec[5][0],0],[vec[6][0],0],[vec[7][0],0],[vec[8][0],0],[vec[9][0],0],
		  [vec[10][0],0],[vec[11][0],0],[vec[12][0],0],[vec[13][0],0],[vec[14][0],0]]
	result = svd.matrixmultiply(mat,matvec)
	return [[result[0][0]],[result[1][0]],[result[2][0]],[result[3][0]],[result[4][0]]]


def f(t,x,u):
    x1,x2,v1,v2,theta = x
    u1,u2 = u
    return np.array([v1,v2,u1*math.cos(theta)-K*v1,u1*math.sin(theta)-K*v2,u2])

def h(t,x,u):
    x1,x2,v1,v2,theta = x
    u1,u2 = u
    return np.array([x1,x2,theta])

def get_u(j):
    u1 = U[2*(4-j)][0]
    u2 = U[2*(4-j)+1][0]
    return np.array([u1,u2])

def sim(f, t, x, u, dt=1e-4):
    j, t_, x_, u_ = 0, [0], [x], [u]
    while j*dt < t:
        t_.append((j+1)*dt)
        u_.append(get_u(int(j/time_for_reaching_subgoal)))
        x1,x2,v1,v2,theta = x_[-1]+dt*f(j*dt,x_[-1],u_[-1])
        theta = theta % (2*math.pi)
        if theta > math.pi:
            theta -= 2*math.pi
        x_.append(np.array([x1,x2,v1,v2,theta]))
        j += 1
    return np.array(t_),np.array(x_)

K = 1
zero = [0,0,0,0,0,0,0,0,0,0]
zero1 = [0,0,0,0,0]
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

print "(1) Observability check:"
CA = svd.matrixmultiply(C,A)
CA2 = svd.matrixmultiply(CA,A)
CA3 = svd.matrixmultiply(CA2,A)
CA4 = svd.matrixmultiply(CA3,A)

Observe_matrix = np.concatenate((CA4,CA3,CA2,CA,C),axis=0) # 15*5
print Observe_matrix
# This is the observability matrix I get:
#[[ 0.          0.         -1.          0.         -0.70710678]
# [ 0.          0.          0.         -1.          0.70710678]
# [ 0.          0.          0.          0.          0.        ]
# [ 0.          0.          1.          0.          0.70710678]
# [ 0.          0.          0.          1.         -0.70710678]
# [ 0.          0.          0.          0.          0.        ]
# [ 0.          0.         -1.          0.         -0.70710678]
# [ 0.          0.          0.         -1.          0.70710678]
# [ 0.          0.          0.          0.          0.        ]
# [ 0.          0.          1.          0.          0.        ]
# [ 0.          0.          0.          1.          0.        ]
# [ 0.          0.          0.          0.          0.        ]
# [ 1.          0.          0.          0.          0.        ]
# [ 0.          1.          0.          0.          0.        ]
# [ 0.          0.          0.          0.          1.        ]]


Observe_rank = matrix_rank(Observe_matrix)
print "rank = " + str(Observe_rank)  # rank = 5 --> fully observable.


print "(2) find steering input sequence"
# we need to transpose first before feeding into svd function
# due to shape constraint
a,b,c = svd.svd(Observe_matrix)
print "a"
print a
print "b"
print b
print "c"
print c
# This is the "a" I get: (15*15)
#
#[[  0.00000000e+00   0.00000000e+00   9.71445147e-17  -3.53553391e-01    -3.77964473e-01]
# [ -0.00000000e+00   0.00000000e+00   1.38777878e-16  -3.53553391e-01     3.77964473e-01]
# [ -0.00000000e+00  -0.00000000e+00   4.29837141e-17  -6.86314164e-18    -8.90930982e-18]
# [ -0.00000000e+00  -0.00000000e+00  -1.24900090e-16   3.53553391e-01     3.77964473e-01]
# [ -0.00000000e+00  -0.00000000e+00  -1.38777878e-16   3.53553391e-01    -3.77964473e-01]
# [ -0.00000000e+00  -0.00000000e+00   0.00000000e+00  -0.00000000e+00    -0.00000000e+00]
# [ -0.00000000e+00  -0.00000000e+00   9.71445147e-17  -3.53553391e-01    -3.77964473e-01]
# [ -0.00000000e+00  -0.00000000e+00   1.38777878e-16  -3.53553391e-01     3.77964473e-01]
# [ -0.00000000e+00  -0.00000000e+00   0.00000000e+00  -0.00000000e+00    -0.00000000e+00]
# [ -0.00000000e+00  -0.00000000e+00   5.00000000e-01   3.53553391e-01     1.88982237e-01]
# [ -0.00000000e+00  -0.00000000e+00  -5.00000000e-01   3.53553391e-01    -1.88982237e-01]
# [ -0.00000000e+00  -0.00000000e+00   0.00000000e+00  -0.00000000e+00    -0.00000000e+00]
# [ -1.00000000e+00  -0.00000000e+00   0.00000000e+00  -0.00000000e+00    -0.00000000e+00]
# [ -0.00000000e+00  -1.00000000e+00   0.00000000e+00  -0.00000000e+00    -0.00000000e+00]
# [ -0.00000000e+00  -0.00000000e+00  -7.07106781e-01  -2.49800181e-16     2.67261242e-01]]


# This is the "b" I get: (15*5)
#
#[1.0, 1.0, 0.9999999999999996, 1.9999999999999996, 2.6457513110645903]


# This is the "c" I get: (5x5)
#
#[[-1.0,  -0.0,                 -0.0,                     0.0,                 0.0], 
# [-0.0,  -1.0,                 -0.0,                     0.0,                 0.0],
# [-0.0,  -0.0,   0.4999999999999998,      0.7071067811865471,  0.5000000000000004],
# [-0.0,  -0.0, -0.50000000000000011,     0.70710678118654768, -0.4999999999999995],
# [-0.0,  -0.0, -0.70710678118654746, -7.0776717819853729e-16, 0.70710678118654768]]

u = np.concatenate((a,[zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero]),axis=1)
w = np.concatenate((np.diag(b),[zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1]),axis=0)
v_star = svd.transpose(c)

print "Now this should be identical to the Observe_matrix:"
print svd.matrixmultiply(u,svd.matrixmultiply(w,v_star)) # Observe_matrix = u * w * v_star



# from the simulation output: (from spaceship_control.py)
# h(t,x,u) = [ 1.01404706  1.01556914  0.94539816]
# h(t,x,u) = [ 1.04199006  1.04898614  0.94539816]
# h(t,x,u) = [ 1.06583733  1.07740244  0.98728607]
# h(t,x,u) = [ 1.08615721  1.10162570  0.94539816]
# h(t,x,u) = [ 1.10345371  1.12225025  0.98728607]

y_bar = [[ 1.01404706],  [1.01556914],  [0.94539816],
         [ 1.04199006],  [1.04898614],  [0.94539816],
         [ 1.06583733],  [1.07740244],  [0.98728607],
         [ 1.08615721],  [1.10162570],  [0.94539816],
         [ 1.10345371],  [1.12225025],  [0.98728607]]

print "y_bar = " + str(y_bar)


v = pinv(v_star)
u_star = pinv(u)
w_dagger = np.diag(np.reciprocal(np.diag(w)))
w_dagger = np.concatenate((w_dagger,[zero,zero,zero,zero,zero]),axis=1)



Xo = matmul(svd.matrixmultiply(svd.matrixmultiply(v,w_dagger),u_star),y_bar)
Xo_real = [[1.],[1.],[0.],[0.],[math.pi/4.]]

print "See the estimation result:"
print "     Xo = " + str(Xo)
print "Real Xo = " + str(Xo_real)
