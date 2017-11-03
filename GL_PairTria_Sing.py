# -*- coding: utf-8 -*-
"""
Created on Sun May 21 16:59:17 2017

@author: adrian
"""


#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 15:23:31 2016

@author: 3417212
"""

import numpy as np
from numpy.linalg import norm
from numpy.linalg import det, inv
import matplotlib.pyplot as plt
import sys



"""
on se place sur [-1;1]x[-1;1]
"""


# regular function
def f(x,y):
    #return x + y    
    return np.exp(- norm(x-y)**2) # exp(x[0]**2 + x[1]**2 + y[0]**2 + y[1]**2) #return sum(exp(x**2 + y**2))  #norm(x)**2 + norm(y)**2	#sin(norm(x) + norm(y))
    
# singular function
def h(x,y):   
    return 1/(norm(x - y))
				

	

	

# for a given triangle we will determine recursively smaller triangles in this triangle
def generate_smaller_triangles(T):
    # T = triangle with points A,B,C

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
    
    
    if depth == 0:
        return [T]
    
    if depth == 1:
        return generate_smaller_triangles(T)
        
    else:

        list_triangles = generate_smaller_triangles(T)
        
        newlist_triangles = []
        for tria in list_triangles:

            # as we know that A = (-1, -1) has a real singularity we continue the maillage only in this direction
            # TO DO: Exclude maybe triangle with the point A
            if any((np.array([-1, -1]) == x).all() for x in tria):
                new_triangles = create_triangles_recursively(depth-1, tria)
                newlist_triangles.extend(new_triangles)
            else:
                newlist_triangles.extend([tria])

        
        return newlist_triangles #list_triangles       
            





def calculateInt(order, t, w, tria1, tria2):
	
     
    a1 = tria1[0]
    b1 = tria1[1]
    c1 = tria1[2]
    
    a2 = tria2[0]
    b2 = tria2[1]
    c2 = tria2[2]

    # calculate integral
    I = 0
    for i, xi in enumerate(t):
        for j, etai in enumerate(t):
            #print "i =", i, "j =", j, "xi =", xi, "etai =", etai
            # transformation - M*vector
            
            x = a1 + (1+xi)/2*(b1-a1) + (1-xi)*(1-etai)/4*(c1-a1)
            y = a2 + (1+xi)/2*(b2-a2) + (1-xi)*(1-etai)/4*(c2-a2)
#            x = a[0] + (1+xi)/2*(b[0]-a[0]) + (1-xi)*(1-etai)/4*(c[0]-a[0])
#            y = a[1] + (1+xi)/2*(b[1]-a[1]) + (1-xi)*(1-etai)/4*(c[1]-a[1])

            #print "x =", x, " y =", y
            
            # regular function
            #I = I + w[i]*w[j]*f(x, y)*(1-xi)/8
			
            # singular function
            I = I + w[i]*w[j]*h(x, y)*(1-xi)/8


    Jacob_tria1 = (b1[0]-a1[0])*(c1[1]-a1[1]) - (c1[0]-a1[0])*(b1[1]-a1[1])
    Jacob_tria2 = (b2[0]-a2[0])*(c2[1]-a2[1]) - (c2[0]-a2[0])*(b2[1]-a2[1])

    #print "Jacobien: ", Jacob_tria
    I = np.abs(Jacob_tria1)*np.abs(Jacob_tria2)*I
    #I = det(M)*I
    #print "I =", I
    return I
 



def calc_Int_Tria(T1, T2, order, depth_tria_rec, t, w):
    
    # calculate recursively all triangles
    listTriangles1 = create_triangles_recursively(depth_tria_rec, T1)
    listTriangles2 = create_triangles_recursively(depth_tria_rec, T2)
    
    
    # calculate for each triangle integral
    sum_int_trias = 0
    i = 0
    for tria1 in listTriangles1:
        for tria2 in listTriangles2:
        #print "tria = ", tria
            
            #CHECK CONDITION
            if (not (any((np.array([-1, -1]) == x).all() for x in tria1))) and (not (any((np.array([-1, -1]) == x).all() for x in tria2))):
                i = i + 1
                #print tria1
                #print tria2
                #print
                int_tria = calculateInt(order, t, w, tria1, tria2)
        
                sum_int_trias = sum_int_trias + int_tria
                


    #print "sum_int_trias =", sum_int_trias
    return sum_int_trias
    
    
  
    
    
def calc_error(T1, T2, order, depth_tria_rec, t, w):
    
    I = calc_Int_Tria(T1, T2, order, depth_tria_rec, t, w)    
    I2 = calc_Int_Tria(T1, T2, order, depth_tria_rec + 1, t, w) 
    print "I =", I
    print "I2 =", I2
    err = abs(I - I2)/abs(I)
    #print "err =", err
    return err






###############################################################################
###############################################################################

###############################################################################
###############################################################################

# in this triangle there is a singularity at the point A = (-1,-1)
A1 = np.array([-1,-1])
B1 = np.array([1,0])
C1 = np.array([-1,2])

# in this triangle there is a singularity at the point A = (-1,-1)
A2 = np.array([-1,-1])
B2 = np.array([-1, -3])
C2 = np.array([-2,-1])

#A1 = np.array([0,0])
#B1 = np.array([1,1])
#C1 = np.array([0,1])
#
#A2 = np.array([1,0])
#B2 = np.array([2,0])
#C2 = np.array([2,1])

#==============================================================================
#A=np.array([0,1])
#B=np.array([0.9,1])
#C=np.array([0,2])
#==============================================================================

T1 = [A1,B1,C1]
T2 = [A2, B2, C2]

#new = create_triangles_recursively(3, T1)
#print(len(new))


order = 4 # ordre for the method gaussLegendre

#numberSubIntervalls = 5

depth_tria_rec = 1 # depth how deep we divise triangles

rangeDepth = range(0, 5)

t,w = np.polynomial.legendre.leggauss(order) # t in [-1,1], w weight, ordre 50


#integr = calc_Int_Tria(T, order, depth_tria_rec, t, w)

#err = calc_error(T1, T2, order, depth_tria_rec, t, w)


#plt.figure()

err = np.array([])
for depth in rangeDepth:
	
     
	err = np.append( err, calc_error(T1, T2, order, depth, t, w))
	
#print "coucou", calculateInt(int1, int2, order, t, w)	
plt.plot(rangeDepth, np.log(err))
plt.title("Singular case")
plt.ylabel("log(error)")
plt.xlabel("depth rec. tria")
#plt.hold(True)
#plt.legend()
plt.show()


