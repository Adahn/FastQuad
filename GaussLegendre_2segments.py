#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 15:23:31 2016

@author: 3417212
"""
import time
import numpy as np
import matplotlib.pyplot as plt
#plt.close("all")

"""
on se place sur [-1;1]x[-1;1]
"""

def somme_de_GaussLegendre2D(N, alpha1, beta1, alpha2, beta2,func ):
    x,y = np.polynomial.legendre.leggauss(N)
    x1 = (beta1 - alpha1)/2 * x + (alpha1 + beta1)/2
    x2 = (beta2 - alpha2)/2 * x + (alpha2 + beta2)/2
    I_exp = (beta1-alpha1)*(beta2-alpha2)*sum(y*func(x1,x2))/4.
    return I_exp
    
def estimateur_erreur2D(N, func, alpha1, beta1, alpha2, beta2 ):
    Sf = somme_de_GaussLegendre2D(N,alpha1, beta1, alpha2, beta2, func)
    #print "Sf = ", Sf, " pour n = ",N
    Sf2 = somme_de_GaussLegendre2D(N, alpha1, (alpha1 + beta1)/2., alpha2, (alpha2 + beta2)/2. , func)  \
        + somme_de_GaussLegendre2D(N, (alpha1 + beta1)/2., beta1, (alpha2 + beta2)/2., beta2, func) \
        + somme_de_GaussLegendre2D(N, (alpha1 + beta1)/2., beta1, alpha2, (alpha2 + beta2)/2., func) \
        + somme_de_GaussLegendre2D(N,  alpha1, (alpha1 + beta1)/2., (alpha2 + beta2)/2., beta2, func)

    #print "Sf2 = ", Sf2, " pour n = ",N
    return (Sf - Sf2)/10.0

def f(x,y):
    return x[1]**2 + x[0]**2 + y[1]**2 + y[0]**2
    
def test(x,y):
    #print("x: ",x, " y: ",y )
    #print(np.array_equal(x,y))
    return y[1]+x[1]+y[0]+x[0]
    
def h(x,y):   
    #print(np.array_equal(x,y))
    #print("norm x : ", np.linalg.norm(x), "norm y : ", np.linalg.norm(y))
    return np.log(sqrt(x[1]**2 + y[1]**2 + (x[0]**2 + y[0]**2) ) )

def segment(a,b,t):
    s1=a[0]+t*(b[0]-a[0])
    s2=a[1]+t*(b[1]-a[1])
    return np.array([s1,s2])

def mesure_intervalle(a,b):
    return np.sqrt( (a[1]-b[1])**2 + (a[0]-b[0])**2 )
    
       
def somme_de_GaussLegendre_segment(N, intervalle1, intervalle2, func ):
    x,y = np.polynomial.legendre.leggauss(N) # x in [-1,1], y constant
    x = ( x+1 )/2. # translate points to [0,1] 
    #print "VALUES:", x[3:5]

    a1=intervalle1[0]
    b1=intervalle1[1]
    a2=intervalle2[0]
    b2=intervalle2[1]
    gamma=np.array([[0,0]])
    sigma=np.array([[0,0]])
    for i in x:
        #print gamma.shape
        #print np.array([segment(a1,b1,i)]).shape
        gamma=np.concatenate((gamma, np.array([segment(a1,b1,i)])), axis = 0)
        sigma=np.concatenate((sigma, np.array([segment(a2,b2,i)])), axis=0)
        #gamma=np.append(gamma, [segment(a1,b1,i)])
        #sigma=np.append(sigma, [segment(a2,b2,i)])
        #print "gamma : ", gamma
        #print "bbb"
#    print("sigma : ", sigma)
    gamma = np.transpose(gamma[1::])
    sigma = np.transpose(sigma[1::])
    #print("sigma : ", sigma)
    #print "shape y ", y.shape
    #print "shape func(gamma,sigma) ", func(gamma,sigma)
    #I_exp=sum(np.dot(y, func(gamma,sigma)))*mesure_intervalle(a1,b1)*mesure_intervalle(a2,b2)/2.
    I_exp = sum( y * func(gamma,sigma))* mesure_intervalle(a1,b1)*mesure_intervalle(a2,b2)/2.
    return I_exp   

def estimateur_erreur_segment(N, intervalle1, intervalle2, func):
    Sf = somme_de_GaussLegendre_segment(N,intervalle1, intervalle2, func)
    print "Sf = ", Sf, " pour n = ",N   
    sup1=fin_intervalle(intervalle1)
    sup2=fin_intervalle(intervalle2)
    inf1=debut_intervalle(intervalle1)
    inf2=debut_intervalle(intervalle2)
    Sf2 = somme_de_GaussLegendre_segment(N,sup1, sup2, func) \
    + somme_de_GaussLegendre_segment(N,inf1, inf2, func) \
    + somme_de_GaussLegendre_segment(N,inf1, sup2, func) \
    + somme_de_GaussLegendre_segment(N,sup1, inf2, func)
    print "Sf2 = ", Sf2, " pour n = ",N
    return abs((Sf - Sf2) / 10.)
    
def debut_intervalle(intervalle):
	temp=[ (x+y)/2. for x,y in zip(intervalle[0],intervalle[1]) ]
	I=[ intervalle[0], temp ]
	return I

def fin_intervalle(intervalle):
	temp=[ (x+y)/2. for x,y in zip(intervalle[0],intervalle[1]) ]
	I=[ temp, intervalle[1] ]
	return I
#

###############################################################################
###############################################################################

epsilon = 40

a1 = np.array([-epsilon, -1])
b1 = np.array([-epsilon, +1])

a2 = np.array([epsilon, -1])
b2 = np.array([epsilon, +1])

intervalle1=[a1, b1]
intervalle2=[a2, b2]

N = range(100,1000,100)
temps=np.array([])
etest=errh=errg = errf = np.array([])

for n in N:
    debut=time.time()
    #etest = np.append(etest, estimateur_erreur_segment(n, intervalle1, intervalle2,test))
    #errf = np.append(errf, estimateur_erreur_segment(n, intervalle1, intervalle2,f))
    errh = np.append(errh, estimateur_erreur_segment(n,intervalle1, intervalle2,h ) )
    fin=time.time()
    tps=fin-debut
    temps=np.append(temps,tps)
    
img=figure(3)
img=subplot(223)
plt.plot(N,temps,"-*")
plt.xlabel("N")
plt.ylabel("temps(sec)")
plt.show()

#print "sigma = ", intervalle1
#print "gamma = " ,intervalle2

#print "de mesure :"
#print mesure_intervalle(intervalle1[0], intervalle1[1])
#print mesure_intervalle(intervalle2[0], intervalle2[1])

#print(etest)
#print(errf)
print(errh)

fig4=plt.figure()
plt.plot(np.log(N), np.log(errh), 'r')
plt.xlabel('log(N)')
plt.ylabel('log(errh)')
fig4.suptitle('erreur de GL pour h(x,y)=ln|x-y| pour epsilon=40', fontsize=20)
fig4.savefig('fGL.jpg')
plt.show()