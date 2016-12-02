from numpy.linalg import matrix_rank
from numpy.linalg import pinv
import numpy as np
import pylab as plt
import math
import svd


def matmul(mat,vec):
	matvec = [[vec[0][0],0],[vec[1][0],0],[vec[2][0],0],[vec[3][0],0],[vec[4][0],0],[vec[5][0],0],[vec[6][0],0],[vec[7][0],0],[vec[8][0],0],[vec[9][0],0],[vec[10][0],0],[vec[11][0],0],[vec[12][0],0],[vec[13][0],0],[vec[14][0],0],[vec[15][0],0],[vec[16][0],0],[vec[17][0],0],[vec[18][0],0],[vec[19][0],0],[vec[20][0],0],[vec[21][0],0],[vec[22][0],0],[vec[23][0],0],[vec[24][0],0],[vec[25][0],0],[vec[26][0],0],[vec[27][0],0],[vec[28][0],0],[vec[29][0],0]]
	result = svd.matrixmultiply(mat,matvec)
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
zero = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
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

# initial condition:  Xo = [1,1,0,0,pi/4]
#                     Uo = [1,1]
# final condition:    Xf = [3,8,0,0,pi]
#                     Uf = [0,0]

print "(1) Observability check:"
CA = svd.matrixmultiply(C,A)
CA2 = svd.matrixmultiply(CA,A)
CA3 = svd.matrixmultiply(CA2,A)
CA4 = svd.matrixmultiply(CA3,A)
CA5 = svd.matrixmultiply(CA4,A)
CA6 = svd.matrixmultiply(CA5,A)
CA7 = svd.matrixmultiply(CA6,A)
CA8 = svd.matrixmultiply(CA7,A)
CA9 = svd.matrixmultiply(CA8,A)

Observe_matrix = np.concatenate((CA9,CA8,CA7,CA6,CA5,CA4,CA3,CA2,CA,C),axis=0) # 30*3
print Observe_matrix
# This is the observability matrix I get:
'''
[[  1.04333183e+13  -2.10375831e+12  -7.59205100e+10]
 [  1.77767749e+12  -3.58447688e+11  -1.29355218e+10]
 [  4.93069031e+13  -9.94216569e+12  -3.58792875e+11]
 [  3.76155896e+11  -7.58474686e+10  -2.73718213e+09]
 [  6.40910742e+10  -1.29230767e+10  -4.66347100e+08]
 [  1.77767749e+12  -3.58447688e+11  -1.29355218e+10]
 [  1.35616660e+10  -2.73454874e+09  -9.86828818e+07]
 [  2.31065801e+09  -4.65898413e+08  -1.68058951e+07]
 [  6.40910742e+10  -1.29230767e+10  -4.66347100e+08]
 [  4.88942008e+08  -9.85879396e+07  -3.55763455e+06]
 [  8.33019493e+07  -1.67897188e+07  -6.05124317e+05]
 [  2.31065801e+09  -4.65898413e+08  -1.68058951e+07]
 [  1.76275788e+07  -3.55421158e+06  -1.28185165e+05]
 [  3.00145169e+06  -6.04541188e+05  -2.14511411e+04]
 [  8.33019493e+07  -1.67897188e+07  -6.05124317e+05]
 [  6.35484503e+05  -1.28061759e+05  -4.61410132e+03]
 [  1.08038936e+05  -2.14301332e+04  -7.49872716e+02]
 [  3.00145169e+06  -6.04541188e+05  -2.14511411e+04]
 [  2.28914803e+04  -4.60965296e+03  -1.62419085e+02]
 [  3.80294980e+03  -7.49118823e+02  -8.37050066e+00]
 [  1.08038936e+05  -2.14301332e+04  -7.49872716e+02]
 [  8.23715160e+02  -1.62258869e+02  -5.65924002e+00]
 [  1.32506500e+02  -8.34400000e+00  -8.00000000e-04]
 [  3.80294980e+03  -7.49118823e+02  -8.37050066e+00]
 [  2.87009000e+01  -5.65350000e+00  -2.00000000e-04]
 [  0.00000000e+00   0.00000000e+00   1.00000000e+00]
 [  1.32506500e+02  -8.34400000e+00  -8.00000000e-04]
 [  1.00000000e+00   0.00000000e+00   0.00000000e+00]
 [  0.00000000e+00   1.00000000e+00   0.00000000e+00]
 [  0.00000000e+00   0.00000000e+00   1.00000000e+00]]
'''


Observe_rank = matrix_rank(Observe_matrix)
print "rank = " + str(Observe_rank)  # rank = 3 --> fully observable.


print "(2) find steering input sequence"
# we need to transpose first before feeding into svd function
# due to shape constraint
a,b,c = svd.svd(Observe_matrix)


# This is the "a" I get: (30*30)
'''
[[ -2.06752633e-01  -6.60343207e-01  -6.34190843e-01]
 [ -3.52274770e-02   4.94296295e-01  -5.02801802e-01]
 [ -9.77093922e-01   1.03529721e-01   1.69907217e-01]
 [ -7.45412132e-03   1.19624786e-02   1.22650501e-02]
 [ -1.27006499e-03   1.77514504e-01   1.74459959e-01]
 [ -3.52274770e-02   4.94296295e-01  -5.02801802e-01]
 [ -2.68745749e-04   5.28837683e-03  -5.77741481e-03]
 [ -4.57892527e-05   2.47250717e-02  -3.15056066e-02]
 [ -1.27006499e-03   1.77514504e-01   1.74459959e-01]
 [ -9.68914976e-06   1.88272669e-03   2.02969877e-03]
 [ -1.65073168e-06   8.62167304e-03   1.13271243e-02]
 [ -4.57892527e-05   2.47250717e-02  -3.15056066e-02]
 [ -3.49317505e-07   2.43499460e-04  -3.82659211e-04]
 [ -5.94758773e-08   9.02839412e-04  -2.30204793e-03]
 [ -1.65073168e-06   8.62167304e-03   1.13271243e-02]
 [ -1.25928316e-08   9.16587596e-05   1.31903994e-04]
 [ -2.13959487e-09   4.22217894e-04   7.37398594e-04]
 [ -5.94758773e-08   9.02839412e-04  -2.30204793e-03]
 [ -4.53606543e-10   8.51691929e-06  -2.76920702e-05]
 [ -7.52907976e-11   2.70098250e-05  -1.63940711e-04]
 [ -2.13959487e-09   4.22217894e-04   7.37398594e-04]
 [ -1.63084469e-11   4.54437878e-06   8.63243253e-06]
 [ -2.55514342e-12   2.15628265e-05   4.87279457e-05]
 [ -7.52907976e-11   2.70098250e-05  -1.63940711e-04]
 [ -5.68211262e-13   2.25339692e-07  -1.95965503e-06]
 [  1.38558935e-16   3.39393405e-07  -1.14771454e-05]
 [ -2.55514342e-12   2.15628265e-05   4.87279457e-05]
 [ -1.90413809e-14   2.35509285e-07   5.72559211e-07]
 [  3.83947385e-15   1.15573023e-06   3.25372161e-06]
 [  1.38558935e-16   3.39393405e-07  -1.14771454e-05]]
'''

# This is the "b" I get: (30*3)
#
#[51479756850486.95, 814769.7748875112, 83729.83925419993]


# This is the "c" I get: (3x3)
#
#[[-0.9802456581022319, 0.19188584704972583, 0.047940290713881256], 
# [0.19765518046376787, 0.94165406071451385, 0.27243358745885632], 
# [0.0071329802717776131, 0.27652748804029625, -0.96097953617679444]]


u = np.concatenate((a,[zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero]),axis=1)
w = np.concatenate((np.diag(b),[zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1,zero1]),axis=0)
v_star = svd.transpose(c)

print "Now this should be identical to the Observe_matrix:"
print svd.matrixmultiply(u,svd.matrixmultiply(w,v_star)) # Observe_matrix = u * w * v_star



# from the simulation output: (from rip_control.py)
'''
h(t,x,u) = [ 0.02387942  0.08783102  0.11384511]
h(t,x,u) = [ 0.04674999  0.08953041  0.22721648]
h(t,x,u) = [ 0.06897661  0.09236759  0.34166801]
h(t,x,u) = [ 0.09087787  0.09636078  0.45869429]
h(t,x,u) = [ 0.11747108  0.10157773  0.58672335]
h(t,x,u) = [ 0.1394019   0.10806533  0.71314037]
h(t,x,u) = [ 0.16193178  0.1158517   0.8468102 ]
h(t,x,u) = [ 0.1853321   0.12501718  0.98931885]
h(t,x,u) = [ 0.20986185  0.13565823  1.14230714]
h(t,x,u) = [ 0.2357641   0.14788802  1.30747784]
'''


y_bar = [[ 0.02387942],[  0.08783102],[  0.11384511],
         [ 0.04674999],[ 0.08953041 ],[ 0.22721648],
         [ 0.06897661],[ 0.09236759],[ 0.34166801],
         [ 0.09087787],[ 0.09636078],[ 0.45869429],
         [ 0.11747108],[ 0.10157773],[ 0.58672335],
         [ 0.1394019],[ 0.10806533],[ 0.71314037],
         [ 0.16193178],[ 0.1158517],[ 0.8468102 ],
         [ 0.1853321],[ 0.12501718],[ 0.98931885],
         [ 0.20986185],[ 0.13565823],[ 1.14230714],
         [ 0.2357641],[ 0.14788802],[ 1.30747784]]

print "y_bar = " + str(y_bar) # 30*1


v = pinv(v_star) # 3*3
u_star = pinv(u) # 30*30
w_dagger = np.diag(np.reciprocal(np.diag(w))) # 3*30
w_dagger = np.concatenate((w_dagger,[zero,zero,zero]),axis=1)

Xo = matmul(svd.matrixmultiply(svd.matrixmultiply(v,w_dagger),u_star),y_bar)
Xo_real = [[0.],[5.*math.pi/180],[0.]]

print "See the estimation result:"
print "     Xo = " + str(Xo)
print "Real Xo = " + str(Xo_real)
