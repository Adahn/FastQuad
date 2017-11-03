#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 15:23:31 2016

@author: 3417212
"""

import numpy as np
import matplotlib.pyplot as plt
#plt.close("all")
import time
"""
on se place sur [-1;1]x[-1;1]
"""
def somme_de_GaussLegendre2D(N, alpha1, beta1, alpha2, beta2,func ):
    x,y = np.polynomial.legendre.leggauss(N)
    
    #homothétie de [-1;1] vers [alpha;beta]
    x1 = (beta1 - alpha1)/2 * x + (alpha1 + beta1)/2
    x2 = (beta2 - alpha2)/2 * x + (alpha2 + beta2)/2
    
    I_exp = (beta1-alpha1)*(beta2-alpha2)*sum(y*func(x1,x2))/2.
#    I_exp = 0
#    for petit_x1 in x1:
#        for petit_x2 in x2:
#            I_exp += y*func(petit_x1, petit_x2)
#    I_exp = I_exp *(beta1-alpha1) *(beta2-alpha2)
    return I_exp
    
def estimateur_erreur2D(N, func, alpha1, beta1, alpha2, beta2 ):
    Sf = somme_de_GaussLegendre2D(N,alpha1, beta1, alpha2, beta2, func)
    print "Sf = ", Sf, " pour n = ",N
    Sf2 = somme_de_GaussLegendre2D(N, alpha1, (alpha1 + beta1)/2., alpha2, (alpha2 + beta2)/2. , func)  \
        + somme_de_GaussLegendre2D(N, (alpha1 + beta1)/2., beta1, (alpha2 + beta2)/2., beta2, func) \
        + somme_de_GaussLegendre2D(N, (alpha1 + beta1)/2., beta1, alpha2, (alpha2 + beta2)/2., func) \
        + somme_de_GaussLegendre2D(N,  alpha1, (alpha1 + beta1)/2., (alpha2 + beta2)/2., beta2, func)

    print "Sf2 = ", Sf2, " pour n = ",N
    return Sf - Sf2

def f(x,y):
    return x**2 + y**2
    
def h(x,y):
    I = np.ones(size(x))
    epsilon =  1e-10
    return np.log( abs(x-y) + epsilon )

###############################################################################
###############################################################################


N = np.linspace(100, 1000, 10)
valh = errh=errg = errf = np.array([])
temps=np.array([])
for n in N:
     debut=time.time()
   # errf = np.append(errf, estimateur_erreur2D(n, f, -1,1,-1,1))
     errh = np.append(errh, estimateur_erreur2D(n, h, -1,1,-1,1))
     fin=time.time()
     tps=fin-debut
     temps=np.append(temps,tps)

"""
fig4=plt.figure()
plt.plot(np.log(N), np.log(errf), 'r')
plt.xlabel('log(N)')
plt.ylabel('log(errf)')
fig4.suptitle('erreur de Gauss Legendre pour f(x,y)=x^2+y^2', fontsize=20)
fig4.savefig('f2D.jpg')
plt.show()
"""
fig3=plt.figure(3)
plt.subplot(221)
plt.title('Rectangle')
plt.plot(N,temps,"-*",label='Gauss-Lengendre')
plt.show()

fig5=plt.figure()
plt.plot(np.log(N), np.log(errh), 'b')
plt.xlabel('log(N)')
plt.ylabel('log(errh)')
plt.show()


M = np.mean(np.log(errh))
print('M =')
print(M)
