from numpy.linalg import matrix_rank
from numpy.linalg import pinv
import numpy as np
import pylab as plt
import math
import svd


def matmul(mat,vec,lenth=3):
	matvec = [[vec[0][0],0],[vec[1][0],0],[vec[2][0],0]]
	result = svd.matrixmultiply(mat,matvec)
	if lenth == 20:
		return [[result[0][0]],[result[1][0]],[result[2][0]],[result[3][0]],[result[4][0]],[result[5][0]],[result[6][0]],[result[7][0]],[result[8][0]],[result[9][0]],[result[10][0]],[result[11][0]],[result[12][0]],[result[13][0]],[result[14][0]],[result[15][0]],[result[16][0]],[result[17][0]],[result[18][0]],[result[19][0]]]
	return [[result[0][0]],[result[1][0]],[result[2][0]]]

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

def get_u(j):
    scale = 1
    u1 = U[2*(4-j)][0]*scale
    u2 = U[2*(4-j)+1][0]*scale
    return np.array([u1,u2])

def sim(f, t, x, u, dt=1e-4):
    j, t_, x_, u_ = 0, [0], [x], [u]
    while j*dt < t:
        t_.append((j+1)*dt)
	u_now = get_u(int(j/time_for_reaching_subgoal))
        u_.append(u_now)
        x1,x2,x3 = x_[-1]+dt*f(j*dt,x_[-1],u_[-1])
	'''
        x2 = x2 % (2*math.pi)
        if x2 > math.pi:
            x2 -= 2*math.pi
	'''
        x_.append(np.array([x1,x2,x3]))
        j += 1
	if j % time_for_reaching_subgoal == 0:
	    h(t,[x1,x2,x3],[u_now[0],u_now[1]])
    return np.array(t_),np.array(x_)


time_for_reaching_subgoal = 100
K = 1
zero = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
zero1 = [0,0,0]
 
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

# initial condition:  Xo = [0,5*pi/180,0]
#                     Uo = [0,0]
# final condition:    Xf = [0,0,0]
#                     Uf = [0,0]

print "(1) Controllability check:"
AB = svd.matrixmultiply(A,B)
A2B = svd.matrixmultiply(A,AB)
A3B = svd.matrixmultiply(A,A2B)
A4B = svd.matrixmultiply(A,A3B)
A5B = svd.matrixmultiply(A,A4B)
A6B = svd.matrixmultiply(A,A5B)
A7B = svd.matrixmultiply(A,A6B)
A8B = svd.matrixmultiply(A,A7B)
A9B = svd.matrixmultiply(A,A8B)

# try 10 steps to reach the goal:
Control_matrix = np.concatenate((A9B,A8B,A7B,A6B,A5B,A4B,A3B,A2B,AB,B),axis=1)
print Control_matrix

# This is the controllability matrix I get: ( 20 x 3 )
#[[  1.16701537e+15  -7.59205100e+10   4.20747941e+13  -2.73718213e+09
#    1.51693596e+12  -9.86828818e+07   5.46904912e+10  -3.55763455e+06
#    1.97174136e+09  -1.28185165e+05   7.10836030e+07  -4.61410132e+03
#    2.56121254e+06  -1.62419085e+02   9.21922439e+04  -5.65924002e+00
#    3.24514869e+03  -2.00000000e-04   1.13069000e+02   0.00000000e+00]
# [  1.98841556e+14  -1.29355218e+10   7.16889036e+12  -4.66347100e+08
#    2.58459248e+11  -1.68058951e+07   9.31788587e+09  -6.05124317e+05
#    3.35791407e+08  -2.14511411e+04   1.20907168e+07  -7.49872716e+02
#    4.28598874e+05  -8.37050066e+00   1.49822439e+04  -8.00000000e-04
#    1.66878500e+02   1.00000000e+00   0.00000000e+00   0.00000000e+00]
# [  5.51520741e+15  -3.58792875e+11   1.98841556e+14  -1.29355218e+10
#    7.16889036e+12  -4.66347100e+08   2.58459248e+11  -1.68058951e+07
#    9.31788587e+09  -6.05124317e+05   3.35791407e+08  -2.14511411e+04
#    1.20907168e+07  -7.49872716e+02   4.28598874e+05  -8.37050066e+00
#    1.49822439e+04  -8.00000000e-04   1.66878500e+02   1.00000000e+00]]

Control_rank = matrix_rank(Control_matrix)
print "rank = " + str(Control_rank)  # rank = 3 --> fully controllable.


print "(2) find steering input sequence"
# we need to transpose first before feeding into svd function
# due to shape constraint
a,b,c = svd.svd(svd.transpose(Control_matrix))

# This is the "a" I get:
#
#[ [-0.9993498654794869, -0.0071272637174278411, -0.035287233572627599], 
#  [6.5012899712774328e-05, 0.048464528315843115, -0.011676998576551103],
#  [-0.036029884101292584, 0.23268457052307792, 0.97024383110883983],
#  [2.3439041678786934e-06, -0.0014493613989748346, -0.0025442942071282142],
#  [-0.0012989955569460915, -0.96927565896308876, 0.23352928298152278],
#  [8.4501712811427175e-08, 0.0027758413655928028, -0.00052072652414927867],
#  [-4.6832568316426981e-05, 0.028993240612419167, 0.050882432811185693],
#  [3.0452613996451146e-09, -0.00022217965847433794, -0.0001118380299298432],
#  [-1.6883932942690281e-06, -0.055515122718841146, 0.010413768036976991],
#  [1.0965078876376009e-10, 0.00016280137990140909, -2.2500039011460544e-05],
#  [-6.0846067370613779e-08, 0.0044438745940387813, 0.0022365997014579613],
#  [3.8895284833134792e-12, -2.0624958601381759e-05, -4.9077597279398991e-06],
#  [-2.1908831435592958e-09, -0.0032559480129530586, 0.00044996725131403058],
#  [1.3589690967859686e-13, 9.9417464119722777e-06, -9.6723686925493854e-07],
#  [-7.7713748404445089e-11, 0.00041251037193751734, 9.8148179023150697e-05],
#  [1.6573543939706876e-15, -1.6437492305340468e-06, -2.1638408348388538e-07],
#  [-2.7151816628455813e-12, -0.00019883118213644572, 1.9343279854198051e-05],
#  [-6.0991822902836171e-18, 6.2725077296475227e-07, -4.1294925793396564e-08],
#  [-3.3050675846831851e-14, 3.287539780542191e-05, 4.327376259434e-06],
#  [-1.7321807743339877e-16, -1.2234691199719545e-07, -9.6111903148014759e-09]]


# This is the "b" I get:
#
#[5644500796982503.0, 1259325.1360411737, 14826343.873286325]


# This is the "c" I get:
#
# [[ -0.20688713165362999,  0.5935538759810035,   0.7777477168491473],
#  [-0.035250400385678095, 0.78990811613768863, -0.61221121954147295],
#  [ -0.97772957612459632, -0.1540745415950856, -0.14249881263884728]]



a = np.concatenate((a,[zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero]),axis=1)

b = np.concatenate((np.diag(b),[zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1]),axis=0)

v_star = svd.transpose(a)
w = svd.transpose(b)
u = c

print "Now this should be identical to the Control_matrix:"
print svd.matrixmultiply(u,svd.matrixmultiply(w,v_star)) # Control_matrix = u * w * v_star

Xo = [[0.],[5.*math.pi/180],[0.]]
Xf = [[0.],[0.],[0.]]

A2 = svd.matrixmultiply(A,A)
A4 = svd.matrixmultiply(A2,A2)
A8 = svd.matrixmultiply(A4,A4)
A10Xo = matmul(svd.matrixmultiply(A8,A2),Xo)


Z = [[Xf[0][0]-A10Xo[0][0]],[Xf[1][0]-A10Xo[1][0]],[Xf[2][0]-A10Xo[2][0]]]

print "Xf = " + str(Xf)
print "A10Xo = " + str(A10Xo)
print "Xf - A10Xo = " + str(Z)

v = pinv(v_star)
u_star = pinv(u)
w_dagger = np.diag(np.reciprocal(np.diag(w)))
w_dagger = np.concatenate((w_dagger,[zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1]),axis=0)

U = matmul(svd.matrixmultiply(svd.matrixmultiply(v,w_dagger),u_star),Z,20)

print "U = " + str(U)


print "See the simulation result:"
plt.figure()
total_t = 10*1e-4*time_for_reaching_subgoal
x0 = [Xo[0][0],Xo[1][0],Xo[2][0]]
u0 = [0,0]
t_,x_ = sim(f,total_t,x0,u0,dt=1e-4)

lineObjects = plt.plot(t_,x_,'-')
plt.xlabel('time (sec)')
plt.ylabel('state')
plt.legend(iter(lineObjects),('x1','x2','x3'),loc=2)

plt.show()

#------------------------------------------------------------
