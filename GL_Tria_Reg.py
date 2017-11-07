#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 15:23:31 2016

@author: Adrian Ahne

Calculates integral over a triangle for a regular function
"""

import numpy as np
from numpy.linalg import norm
from numpy.linalg import det, inv
import matplotlib.pyplot as plt
import sys
import time


# regular function
def f(x,y):
    #return x + y    
    return np.exp(- norm(x-y)**2) # exp(x[0]**2 + x[1]**2 + y[0]**2 + y[1]**2) #return sum(exp(x**2 + y**2))  #norm(x)**2 + norm(y)**2	#sin(norm(x) + norm(y))
    
# singular function
def h(x,y):   
    return 1/(norm(x - y))
				

	

def generate_smaller_triangles(T):
    """ T = triangle with points A,B,C 
        for a given triangle we will determine recursively 
        smaller triangles in this triangle
    """

    A = T[0]
    B = T[1]
    C = T[2]

    halfAB = (A+B)/2.0
    halfAC = (A+C)/2.0
    halfBC = (B+C)/2.0
    
    tria1 = [A, halfAB, halfAC]    
    tria2 = [halfAB, B, halfBC]
    tria3 = [halfAB, halfAC, halfBC]
    tria4 = [halfAC, halfBC, C]
    
    list_triangles = [tria1, tria2, tria3, tria4]

    return list_triangles    


def create_triangles_recursively(depth, T):
    """ create the trinagles in a recursiv manner 
        getting smaller and smaller towards the singularity
    """
    
    if depth == 0:
        return [T]
    
    if depth == 1:
        return generate_smaller_triangles(T)
        
    else:

        list_triangles = generate_smaller_triangles(T)
        
        newlist_triangles = []
        for tria in list_triangles:
            new_triangles = create_triangles_recursively(depth-1, tria)

               
            newlist_triangles.extend(new_triangles)

        
        return newlist_triangles     



def calculateInt(order, t, w, tria, case):
    """ calculate integral over a triangles """    

    a = tria[0]
    b = tria[1]
    c = tria[2]

    # calculate integral
    I = 0
    for i, xi in enumerate(t):
        for j, etai in enumerate(t):

            x = a[0] + (1+xi)/2*(b[0]-a[0]) + (1-xi)*(1-etai)/4*(c[0]-a[0])
            y = a[1] + (1+xi)/2*(b[1]-a[1]) + (1-xi)*(1-etai)/4*(c[1]-a[1])
            
            # regular function
            if case == "regular":
                I = I + w[i]*w[j]*f(x, y)*(1-xi)/8
			
            # singular function
            elif case == "singular":            
                I = I + w[i]*w[j]*h(x, y)*(1-xi)/8

            else:
                print("error: ", case)

    Jacob_tria = (b[0]-a[0])*(c[1]-a[1]) - (c[0]-a[0])*(b[1]-a[1])
    I = np.abs(Jacob_tria)*I

    return I
 



def calc_Int_Tria(T, order, depth_tria_rec, t, w, case):
    """ calculate integral over all triangles """    

    # calculate recursively all triangles
    listTriangles = create_triangles_recursively(depth_tria_rec, T)
    
    # calculate for each triangle integral
    sum_int_trias = 0
    for tria in listTriangles:
        int_tria = calculateInt(order, t, w, tria, case)
        
        sum_int_trias = sum_int_trias + int_tria

    return sum_int_trias
    
    
  
    
    
def calc_error(T, order, depth_tria_rec, t, w, case):
    """ Calculate relativve error of integral approximation 
        by considering two consequent depth for the recursive
        triangle creation
    """

    I = calc_Int_Tria(T, order, depth_tria_rec, t, w, case)    
    I2 = calc_Int_Tria(T, order, depth_tria_rec + 1, t, w, case) 
    print("I =", I)
    print("I2 =", I2)
    err = abs(I - I2)/abs(I)

    return err






###############################################################################
###############################################################################
# TESTING
###############################################################################
###############################################################################

#A = np.array([0,0])
#B = np.array([0.5,0])
#C = np.array([0,0.5])


#==============================================================================
A=np.array([0,1])
B=np.array([0.9,1])
C=np.array([0,2])
#==============================================================================

T = [A,B,C]

order = 4 # ordre for the method gaussLegendre

depth_tria_rec = 1 # depth how deep we divise triangles

rangeDepth = range(0,6)

t,w = np.polynomial.legendre.leggauss(order) # t in [-1,1], w weight, ordre 50


errReg = np.array([])
for depth in rangeDepth:
	
  #debut=time.time()
  errReg = np.append( errReg, calc_error(T, order, depth, t, w, "regular"))
  #fin=time.time()
  #temps=np.append(temps,fin-debut)


errSing = np.array([])
for depth in rangeDepth:
  errSing = np.append( errSing, calc_error(T, order, depth, t, w, "singular"))


plt.figure()

# sp1
plt.subplot(211)
plt.plot(rangeDepth, np.log(errReg))
plt.title("Regular case")
plt.ylabel("log(error)")
plt.xlabel("depth rec tria")


plt.subplots_adjust(hspace = 0.6)
# sp2
plt.subplot(212)
plt.plot(rangeDepth, np.log(errSing))
plt.title("Singular case")
plt.ylabel("log(error)")
plt.xlabel("depth rec tria")

plt.show()


#print calculateInt(int1, int2, o rder, t, w)	
#plt.figure()
#.plot(rangeDepth, np.log(errReg))
#plt.title("Regular case")
#plt.ylabel("log(error)")
#plt.xlabel("log(depth rec tria)")
#plt.hold(True)
#plt.legend()
#plt.show()

#image=figure(3)
#plt.subplot(224)
#plt.plot(rangeDepth,temps,"b--",linewidth=5)
#plt.show()