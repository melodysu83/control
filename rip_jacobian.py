import sympy as sp
from sympy import *
import numpy as np
import math

X1_value = 0
X2_value = 0
X3_value = 0
U1_value = 0
U2_value = 0

x1,x2,x3,u1,u2 = sp.symbols('x1,x2,x3,u1,u2', real = True)
J = Function('J')

f1 = -(3.8462*x1+0.2166*x3**2*sin(x2)-19.5258*cos(x2)*sin(x2)+0.00011352*x3*cos(x2)+0.1083*x1*x3*sin(x2)-0.1598*x1**2*cos(x2)**2*sin(x2)+0.1083*x1*x3*cos(x2)*sin(x2))/(0.1083*sin(x2)-0.3197*cos(x2)**2+1)
f2 = x3                               
f3 = -(0.000524*x3-90.1468*sin(x2)-9.7629*sin(x2)**2+5.6766*x1*cos(x2)+0.00005676*x3*sin(x2)-0.0799*x1**2*cos(x2)*sin(x2)**2-0.7379*x1**2*cos(x2)*sin(x2)+0.3197*x3**2*cos(x2)*sin(x2)+0.0160*x1*x3*cos(x2)*sin(x2)+0.1598*x1*x3*cos(x2)**2*sin(x2))/(0.1083*sin(x2)-0.3197*cos(x2)**2+1)

g1 = 76.9231/(0.1083*sin(x2)-0.3197*cos(x2)**2+1)
g2 = 0
g3 = 113.5308*cos(x2)/(0.1083*sin(x2)-0.3197*cos(x2)**2+1)

F1 = f1+g1*u1
F2 = f2+g2*u1
F3 = f3+g3*u1

F = sp.Matrix([F1,F2,F3])
A = F.jacobian([x1,x2,x3])
A_sub = F.jacobian([x1,x2,x3]).subs([(x1,X1_value), (x2,X2_value), (x3,X3_value), (u1,U1_value), (u2,U2_value)])

B = F.jacobian([u1,u2])
B_sub = F.jacobian([u1,u2]).subs([(x1,X1_value), (x2,X2_value), (x3,X3_value), (u1,U1_value), (u2,U2_value)])


C = [[1,0,0],[0,1,0],[0,0,1]]
C_sub = C
D = [[0,0],[0,0],[0,0]]
D_sub = D

print "\nLinearized at state Xo: [x1,x2,x3] = " + str([X1_value,X2_value,X3_value])
print "              input Uo: [u1,u2] = " + str([U1_value,U2_value])

print "\nmatrix A"
print A
print "\nmatrix A substituting with initial condition (Xo,Uo)"
print A_sub
print "\nmatrix B"
print B
print "\nmatrix B substituting with initial condition (Xo,Uo)"
print B_sub
print "\n\nmatrix C"
print C
print "\nmatrix C substituting with initial condition (Xo,Uo)"
print C_sub
print "\nmatrix D"
print D
print "\nmatrix D substituting with initial condition (Xo,Uo)"
print D_sub
