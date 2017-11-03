# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 15:27:20 2016

@author: 3100965
"""
import random 
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.linalg import norm

N1=1000 #borne inf
N = 2000 #borne sup
pas=10 #pas
#on choisi aléatoirement des points sur chaque segment
def monteCarlo2Seg(a, b, N):
    t=np.random.uniform(0,1,N)
    
    # parametrisation
    x1 = (1-t)*a[0]+t*b[0]
    x2 = (1-t)*a[1]+t*b[1]        

    point = np.array([x1, x2])                    
   
    return point

#regular function that we want to integrate 
def regular_function(Point1,Point2):
    x1, y1 = Point1
    x2, y2 = Point2
    f=np.exp(-((x1-x2)**2+(y1-y2)**2)) # e^(-|x-y|²)
    return f

#singular function that we want to integrate
def singular_function(Point1,Point2):
    x1, y1 = Point1
    x2, y2 = Point2
    f=1/np.sqrt((x1-x2)**2+(y1-y2)**2) # 1/(|x-y|)
    return f
    
#on choisi aléatoirement des points sur chaque segment
def monteCarlo2Seg(a, b, N):
    t=np.random.uniform(0,1,N)

    # parametrisation
    x1 = (1-t)*a[0]+t*b[0]
    y2 = (1-t)*a[1]+t*b[1]        
    point = np.array([x1, y2])                    

    return point

#on mesure la taille des segments

def mesure_x_y(a1,b1,a2,b2):
    V=norm(b1-a1)*norm(b2-a2)
    return V


#regular : e^(-|x-y|²)
def monteCarlo_err_reg(Point1, Point2):
    x1, y1 = Point1
    x2, y2 = Point2
    N = len(x1)
    
    f_moy = sum(regular_function(Point1,Point2))/N
    f2=sum(regular_function(Point1,Point2)**2)
    Var = np.abs(f2 - N*f_moy**2)/(N-1)
    err = np.sqrt(Var/N)*1.96/f_moy
    
    return err
    
def Int_reg(Point1,Point2,V):
    x1, y1 = Point1
    x2, y2 = Point2
    N = len(x1)
    
    f_moy = sum(regular_function(Point1,Point2))/N
    I=f_moy*V
    return I


# singular : 1/(|x-y|)
def monteCarlo_err_sing(Point1, Point2):
    x1, y1 = Point1
    x2, y2 = Point2
    N=len(x1)
    
    f_moy = sum(singular_function(Point1,Point2))/N
    f2=sum(singular_function(Point1,Point2)**2)
    Var = np.abs((f2-N*f_moy**2)/(N-1))
    err = np.abs(np.sqrt(Var/N)*1.96/f_moy)
    return err

def Int_sing(Point1,Point2,V):
    x1, y1 = Point1
    x2, y2 = Point2
    N = len(x1)
    
    f_moy = sum(singular_function(Point1,Point2))/N
    I=V*f_moy
    return I
    
listEps = np.array([0,0.5, 1, 1.5, 2])

print("===== Case singular =====")
figure()
for eps in listEps:

    a1 = np.array([-eps , -1])
    b1 = np.array([-eps, 1])
    a2 = np.array([eps,- 1])
    b2 = np.array([eps, 1])
    
    err=np.array([])
    for n in range(N1,N,pas):
        err_=0
        for k in range(10):
            Point1 = monteCarlo2Seg(a1, b1, n) # points of segment 1
            Point2 = monteCarlo2Seg(a2, b2, n) # points of segment 2
            err_=err_+monteCarlo_err_sing(Point1,Point2)
        err_=err_/10
        err=np.append(err,err_)
    
    plt.plot(log(range(N1,N,pas)),log(err), label = str(eps))
    plt.plot(log(range(N1,N,pas)), -0.5*log(range(N1,N,pas))+1.55, "k")
    print "eps :",eps
    print "err",err1[10]

plt.legend(loc = 2)
plt.title("Case 1 Segments Singular case $-ln(|x-y|)$ more eps|| pente = -0.5")
plt.xlabel("log(N)")
plt.ylabel("log(err)")

print("===== Case regular =====")
figure()
for eps in listEps:
    err1=np.array([])
    for n in range(N1,N,pas):
        err_=0
        for k in range(10):
            Point1 = monteCarlo2Seg(a1, b1, n) # points of segment 1
            Point2 = monteCarlo2Seg(a2, b2, n) # points of segment 2
            err_=err_+monteCarlo_err_reg(Point1,Point2)
        err_=err_/10
        err1=np.append(err1,err_)
    
    print "eps :",eps
    print "err",err1[10]
    plt.plot(log(range(N1,N,pas)),log(err1), label = str(eps))
    #plt.plot(log(range(100,N,10)), -0.5*log(range(100,N,10))-0.55, "r")


plt.legend(loc = 2)
plt.title("Case 2 Segments Regular case $exp(-|x-y|^2)$ more eps|| pente = -0.5")
plt.xlabel("log(N)")
plt.ylabel("log(err)")
