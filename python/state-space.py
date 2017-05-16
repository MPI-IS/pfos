#!/usr/bin/env python
# -*- coding: utf-8

# The code in this file is the computer algebra implementation
# of the probabilistic Nordsieck method. It demonstrates the
# equivalence of the IWP(q) state-space filter and q-step
# q-order RK code.
#
# Schober, M., Särkkä, S., and Hennig, P.: A Probabilistic
# Model for the Numerical Solution of Initial Value Problems
#
# Date: 2016-10-14


## Setup, Configuration

import sympy as sp
sp.init_printing(use_unicode=True)

## Helper functions

# creates a 1x1 matrix with element sym for compatibility
def mat(sym):
	return sp.Matrix([[sym]])

# creates e_i in R^(q+1)
def basisVector(i,q):
	return sp.Matrix(q+1, 1, 
		             lambda m,n: 1 if m == i else 0)

# creates a vector in R^(q+1) with symbolic entries with
# variable name 'ident'
def varVector(q,ident):
	return sp.Matrix(q+1, 1,
		             lambda i,j: sp.symbols(ident + '_%d' % i))


# creates a matrix in R^(q+1)x(q+1) with symbolic entries
# with variable name 'ident'
def varMatrix(q,ident):
	return sp.Matrix(q+1, q+1,
		             lambda i,j: sp.symbols(ident + '_%d%d' % (i,j)))

# creates a symmetric matrix in R^(q+1)x(q+1) with symbolic
# entries with variable name 'ident'
def symVarMatrix(q,ident):
	return sp.Matrix(q+1, q+1,
		             lambda i,j: sp.symbols(ident + '_%d%d' % 
		             	                    (min(i,j), max(i,j))))

# applies limit elementwise
def limitElements(A, t, lim):
	return sp.Matrix(A.rows, A.cols,
		             lambda i,j: sp.limit(A[i,j], t, lim))

# displays a variable during execution
def show(var):
	sp.pprint(var)
	print

## Definitions

h, a, u, v = sp.symbols('h alpha u v')
t0, y0 = sp.symbols('t_0 y_0')

tau, s = sp.symbols('tau sigma')

z = sp.Matrix(1, 10, lambda i,j: sp.symbols('z_%d' % j))

# creates the feedback matrix F for the IWP(q)
def feedbackMatrix(q):
	return sp.Matrix(q+1, q+1, 
		             lambda i,j: 1 if j == i+1 else 0)

# creates the diffusion matrix L for the IWP(q)
def diffusionMatrix(q):
	return sp.Matrix(q+1, 1,
		             lambda i,j: 1 if i == q else 0)

# creates the rescaling matrix B for the IWP(q) and
# step size h
def rescalingMatrix(q,h = h):
	return sp.Matrix(q+1, q+1,
		             lambda i,j: h**i/sp.factorial(i)
		                         if i == j else 0)

# creates the matrix [F, L*L.T; 0, -F.T], needed for
# the discretization
def Phi(F,L):
	return F.row_join(L*L.T).col_join(
		       (0 * F).row_join(-F.T))

# creates the discrete A(h) and Q(h)
def discreteDiffusion(F, L, h = h, ssq = s**2):
	q = F.rows - 1
	P = sp.exp(Phi(F,L)*h)
	A = P[:q+1,:q+1]
	Q = ssq * P[:q+1,q+1:] * A.T
	return (A, Q)

## equivalence of IWP(2) and second-order RK methods
q  = 2

F  = feedbackMatrix(q)
L  = diffusionMatrix(q)
H0 = basisVector(0,q).T
H1 = basisVector(1,q).T

hns = [0, h*a, h*(1-a)]

ms = []
Cs = []

# initialization
mm1p      = sp.zeros(q+1, 1)
(_, Cm1p) = discreteDiffusion(F, L, tau)

# add initial value and first derivative
Km1 = Cm1p * H0.T * (H0 * Cm1p * H0.T)**-1
m   = mm1p + Km1 * (mat(y0) - H0 * mm1p)
C   = (sp.eye(q+1) - Km1 * H0) * Cm1p

for n in range(len(hns)):
	# Predict
	hn = hns[n]
	(A, Q) = discreteDiffusion(F, L, hn)
	mp = A * m
	Cp = A * C * A.T + Q
	# Evaluate
	mlim = sp.limit(mp[0], tau, sp.oo)
	print "Evaluation at"
	show(sp.collect(sp.simplify(mlim),h))
	# Update
	Kn = Cp * H1.T * (H1 * Cp * H1.T)**-1
	m  = mp + Kn * (mat(z[n]) - H1 * mp)
	C  = (sp.eye(q+1) - Kn * H1) * Cp
	# ms.append(mp)
	# Cs.append(Cp)

## equivalence of IWP(3) and third-order RK methods
q  = 3

F  = feedbackMatrix(q)
L  = diffusionMatrix(q)
H0 = basisVector(0,q).T
H1 = basisVector(1,q).T

hns = [0, h*u, h*(v-u), h*(1-v)]

ms = []
Cs = []

# initialization
mm1p      = sp.zeros(q+1, 1)
(_, Cm1p) = discreteDiffusion(F, L, tau)

# add initial value and first derivative
Km1 = Cm1p * H0.T * (H0 * Cm1p * H0.T)**-1
m   = mm1p + Km1 * (mat(y0) - H0 * mm1p)
C   = (sp.eye(q+1) - Km1 * H0) * Cm1p

for n in range(len(hns)):
	# Predict
	hn = hns[n]
	(A, Q) = discreteDiffusion(F, L, hn)
	mp = A * m
	Cp = A * C * A.T + Q
	# Evaluate
	# print "Limiting ..."
	mlim = sp.limit(sp.simplify(mp[0]), tau, sp.oo)
	# print "Simplifying ..."
	# msim = sp.simplify(mlim)
	print "Evaluation at"
	show(sp.collect(sp.simplify(mlim), h))
	# Update
	Kn = Cp * H1.T * (H1 * Cp * H1.T).inv()
	m  = mp + Kn * (mat(z[n]) - H1 * mp)
	C = (sp.eye(q+1) - Kn * H1) * Cp
	# ms.append(mp)
	# Cs.append(Cp)

