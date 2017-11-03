# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 15:27:20 2016

@author: 3100965
"""
import random 
import numpy as np
import math
import matplotlib.pyplot as plt
import time
from scipy.linalg import norm

N1=100 #borne inf
N = 1000 #borne sup
pas=10 #pas
eps = 0.5

a1 = np.array([-eps , -1])
b1 = np.array([-eps, 1])
a2 = np.array([eps, -1])
b2 = np.array([eps, 1])

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
'''    
# singular : -ln(|x-y|)
def monteCarlo_err_sing(Point1, Point2):
    x1, y1 = Point1
    x2, y2 = Point2
    N=len(x1)
    
    f_moy = sum(-np.log(np.sqrt((x1 - x2)**2 + (y1 - y2)**2)))/N
    f2=sum(-np.log(np.sqrt((x1-x2)**2 + (y1-y2)**2))**2)
    Var = np.abs((f2-N*f_moy**2)/(N-1))
    err = np.abs(np.sqrt(Var/N)*1.96/f_moy)

    return err

'''
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
    
#--------------------------------------------------
tps = np.array([])

V=mesure_x_y(a1,b1,a2,b2)
#cas singulier
I=0
Int_=0
cpt=0
err1 = np.array([])
for n in range(N1,N,pas):
    cpt=cpt+1
    err_=0
    for k in range(10):
        Point1 = monteCarlo2Seg(a1, b1, n) # points of segment 1
        Point2 = monteCarlo2Seg(a2, b2, n) # points of segment 2
        err_=err_+monteCarlo_err_sing(Point1,Point2)
        I=I+Int_sing(Point1,Point2,V)
    I=I/10
    Int_=Int_+I
    err_=err_/10
    err1=np.append(err1,err_)
Int_=Int_/cpt
print "Int Sing",Int_

figure()
plt.plot(np.log(range(N1,N,pas)),np.log(err1),"b")
plt.plot(np.log(range(N1,N,pas)), -0.5*np.log(range(N1,N,pas))-2.25, "r")
plt.title("Case 1 Segments Singular case $-ln(|x-y|)$|| pente = -0.5")
plt.xlabel("log(N)")
plt.ylabel("log(err)")
plt.show()

#cas régulier
I=0
Int_=0
cpt=0
err2 = np.array([])

for n in range(N1,N,pas):  
    debut=time.time()
    cpt=cpt+1
    err_=0
    for k in range(10):
        Point3 = monteCarlo2Seg(a1, b1, n) # points of segment 1
        Point4 = monteCarlo2Seg(a2, b2, n) # points of segment 2
        err_=err_+monteCarlo_err_reg(Point3,Point4)
        I=I+Int_reg(Point3,Point4,V)
    I=I/10
    Int_=Int_+I
    #print "I",I
    err_=err_/10
    err2=np.append(err2,err_)
    fin=time.time()
    #print "time ", fin-debut
    tps = np.append(tps, fin-debut)

Int_=Int_/cpt
print "Int Reg",Int_

figure()
plt.plot(np.log(range(N1,N,pas)),np.log(err2),"b")
plt.plot(np.log(range(N1,N,pas)), -0.5*np.log(range(N1,N,pas))+1.05, "r")
plt.title("Case 2 Segments Regular $exp(-|x-y|^2)$|| pente = -0.5")
plt.xlabel("log(N)")
plt.ylabel("log(err)")
plt.show()

figure(3)
subplot(223)
plt.plot(range(N1,N,pas), tps,"-*r",linewidth=5)
plt.title("2 segments")
plt.xlabel("N")
plt.ylabel("temps(sec)")
plt.show()

