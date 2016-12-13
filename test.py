import numpy as np
import control

A = [[-5.65368219902984, 28.7017492282816, -0.000166867558430104], 
     [                0,                0,                     1], 
     [-8.34425988534470, 132.510363075114, -0.000770248419814788]]

C = [[1.,0.,0.],
     [0.,1.,0.],
     [0.,0.,1.]]

L_T = control.place(np.transpose(A),np.transpose(C),[-12.87243308,-3.64492168,-10.86290231])
L_matrix = np.transpose(L_T)
print "\nL matrix = \n" + str(L_matrix)

A_LC = A - np.dot(L_matrix,C)
print "\nA - LC = \n" + str(A_LC)

eigval_newO,eigvec_newO = np.linalg.eig(A_LC)
print "\neigen values of (A_LC) = " + str(eigval_newO)

'''
A_LC=[[  -6.1912818,    28.39236495,   -3.36097806],
      [  -0.30938428,   -0.17804818,   -0.93412002],
      [ -11.70507107,  130.57624305,  -21.01092709]]
print "\nA - LC = \n" + str(A_LC)
'''

next = A_LC
err = [[0],[0.1],[0]]
for i in range(0,10):
	print "iteration " + str(i) + ":"
	next = np.dot(A_LC,np.dot(A_LC,next))
	eigval_newO,eigvec_newO = np.linalg.eig(next)
	print "eigen values of (A_LC) = " + str(eigval_newO)
	print "e = " + str(np.dot(next,err)) + "\n"
