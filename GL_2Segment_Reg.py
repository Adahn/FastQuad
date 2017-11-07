#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 15:23:31 2016

@author: Adrian Ahne
"""

import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
import sys



# regular function
def f(x,y):
	return np.exp(- norm(x-y)**2) 

# singular function
def h(x,y):   
    return 1/(norm(x - y))
				
				

def segment(a,b,t):
    """ Transform points t on interval [-1, 1]"""
    s1 = (b[0]-a[0])/2*t + (b[0]+a[0])/2
    s2 = (b[1]-a[1])/2*t + (b[1]+a[1])/2
    return np.array([s1,s2])


def length_Segment(a,b):
	return norm(b-a) 
    
       


def calculateInt(int1, int2, order, t, w):
	""" Calculates integral over the two segments int1, int2 
		with Gauss-Legendre method
	""" 

	a1, b1 = int1
	a2, b2 = int2	
	
	# calculate parametrisation of vector t
	gamma = segment(a1, b1, t)
	sigma = segment(a2, b2, t)

	I = 0
	for i, gi in enumerate(gamma.T):
		for j, sj in enumerate(sigma.T):
			
			# regular function
			I = I + w[i]*w[j]*f(gi, sj)
			
			# singular function
			#I = I + w[i]*w[j]*h(gi, sj)

	I = length_Segment(a1, b1)*length_Segment(a2, b2)*I

	return I
	


def calculateCompositeInt(int1, int2, Nb_Int, order, t, w):
	""" Composite Gauss-Legendre methode """	

	a1, b1 = int1
	a2, b2 = int2	
	
	I = 0
	for i in range(0, Nb_Int):
		for j in range(0, Nb_Int):
			newInt1 = [np.array(a1 + i*(b1-a1)/Nb_Int), np.array(a1 + (i+1)*(b1-a1)/Nb_Int)]
			newInt2 = [np.array(a2 + j*(b2-a2)/Nb_Int), np.array(a2 + (j+1)*(b2-a2)/Nb_Int)]

			calcI = calculateInt(newInt1, newInt2, order, t, w)
			I = I + calcI
			
	return I
			

	
	
def calculateError(int1, int2, Nb_Int, order, t, w):
	""" Calculate relative error of integral approximation """

	a1, b1 = int1
	a2, b2 = int2

	
	I = calculateCompositeInt(int1, int2, Nb_Int, order, t, w)
	
	# subdivise both segments to calculate error
	
	# segment [a1, b1] becomes [a1, (a1+b1)/2] and [(a1+b1)/2, b1]
	int11 = [int1[0], (int1[0]+int1[1])/2]
	int12 = [(int1[0]+int1[1])/2, int1[1]]
	
	# segment [a2, b2] becomes [a2, (a2+b2)/2] and [(a2+b2)/2, b2]
	int21 = [int2[0], (int2[0]+int2[1])/2]
	int22 = [(int2[0]+int2[1])/2, int2[1]]
	
	# calculate sum of subIntervalls
	I1 = calculateCompositeInt(int11, int21, Nb_Int, order, t, w)
	I2 = calculateCompositeInt(int11, int22, Nb_Int, order, t, w) 
	I3 = calculateCompositeInt(int12, int21, Nb_Int, order, t, w)
	I4 = calculateCompositeInt(int12, int22, Nb_Int, order, t, w)
	I_sub = I1 + I2 + I3 + I4

	err = abs(I - I_sub)/abs(I)

	return err
	
	





###############################################################################
###############################################################################
# TESTING
###############################################################################
###############################################################################


order = 4 # ordre for the method gaussLegendre

numberSubIntervalls = 25

#N = 300 # number of points
# test if we can calculate a straight number for the interval like 200, 242 and NOT 233,343
#if N%order != 0:
#	print "ERROR: Give N and order s.t. N/order is an integer!"
#	exit()

rangeN = range(1, numberSubIntervalls)



t,w = np.polynomial.legendre.leggauss(order) # t in [-1,1], w weight, ordre 50

plt.figure()
for epsilon in [0.1, 1, 2]:
	
	a1 = np.array([-epsilon, -1.0])
	b1 = np.array([-epsilon, +1.0])

	a2 = np.array([epsilon, -1.0])
	b2 = np.array([epsilon, +1.0])

	int1=[a1, b1]
	int2=[a2, b2]
	
	errSing = np.array([])
	for n in rangeN:
		errSing = np.append( errSing, calculateError(int1, int2, n, order, t, w))
	

	plt.plot(np.log(rangeN), np.log(errSing), label=epsilon )
	plt.hold(True)
	plt.legend()

plt.title("regular function for different eps")
plt.ylabel("log(error)")
plt.xlabel("log(number sub intervalls)")

