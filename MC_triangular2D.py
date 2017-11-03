# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 14:39:42 2016

@author: 3100965
"""
import random 
import numpy as np
from numpy.linalg import norm
import math
import matplotlib.pyplot as plt
import time
N1=100 #borne inf
pas=10 #pas
N=1000
#regular function that we want to integrate 
def regular_function(Point1):
    x1, y1 = Point1
    f=np.exp(-(np.abs(x1-y1)**2)) # e^(-|x-y|²)
    return f

#singular function that we want to integrate
def singular_function(Point1):
    x1, y1 = Point1
    f=1/np.abs(x1-y1) # 1/(|x-y|)
    return f
    
def Tria(A,B,C,N):
    r1 = np.random.uniform(0, 1, N)
    r2 = np.random.uniform(0, 1, N)
    r1min=np.minimum(r1,r2)
    r2max=np.maximum(r1,r2)
    x=r1min*A[0]+(r2max-r1min)*B[0]+(1-r2max)*C[0]
    y=r1min*A[1]+(r2max-r1min)*B[1]+(1-r2max)*C[1]
    P=np.array([x,y])     
    return P
    
#formule de Héron
def Aire_Tria(A,B,C):
    a=np.sqrt((B[0]-C[0])**2+(B[1]-C[1])**2)
    b=np.sqrt((A[0]-C[0])**2+(A[1]-C[1])**2)
    c=np.sqrt((B[0]-A[0])**2+(B[1]-A[1])**2)
    P=a+b+c
    p=P/2
    S=np.sqrt(p*(p-a)*(p-b)*(p-c))
    return S
    
#exp(-|x-y|²)
def monteCarlo_Tria_Reg(P,N):
    f_moy=sum(regular_function(P))/N
    f_moy2=sum(regular_function(P)**2)
    Var=np.abs(f_moy2-N*f_moy**2)/(N-1)
    err = np.sqrt(Var/N)*1.96/f_moy
    return err

def Int_reg(Point1,V):
    x1, y1 = Point1
    N = len(x1)
    f_moy=sum(regular_function(Point1))/N
    I=V*f_moy
    return I

#1/|x-y|
def monteCarlo_Tria_Sing(P,N):
    f_moy=sum(singular_function(P))/N
    f_moy2=sum(singular_function(P)**2)
    Var=(f_moy2-N*f_moy**2)/(N-1)
    err = np.sqrt(Var/N)*1.96/f_moy
   # print "f_err", err   
    return err

def Int_sing(Point1,V):
    x1, y1 = Point1
    N = len(x1)
    f_moy=sum(singular_function(Point1))/N
    I=V*f_moy
    return I
    
    
A=np.array([0,1])
B=np.array([0.9,1])
C=np.array([0,2])
Point=Tria(A,B,C,N)
Aire=Aire_Tria(A,B,C)

#plt.figure
#plt.scatter(A[0], A[1])
#plt.scatter(B[0], B[1])
#plt.scatter(C[0], C[1])
#plt.scatter(Point[0],Point[1])
#plt.show()

#
#N=size(Point[0])
temps=np.array([])
Point=np.array([])
err1=np.array([])
err2=np.array([])
figure()
I_Reg=0
I_Sing=0
Int_R=0
Int_S=0
cpt=0
for n in range(100,1000,10):
    debut=time.time()
    cpt=cpt+1
    err1_=0
    err2_=0
    for k in range(10):
        Point=Tria(A,B,C,n)
        err1_=err1_+monteCarlo_Tria_Reg(Point,n)
        err2_=err2_+monteCarlo_Tria_Sing(Point,n)
        I_Reg=I_Reg+Int_reg(Point,V)
        I_Sing=I_Sing+Int_sing(Point,V)
    I_Reg=I_Reg/10
    I_Sing=I_Sing/10
    Int_R=Int_R+I_Reg
    Int_S=Int_S+I_Sing
    err1_=err1_/10
    err2_=err2_/10
    err1=np.append(err1,err1_)
    err2=np.append(err2,err2_)
    fin=time.time()
    temps=np.append(temps,fin-debut)
Int_R=Int_R/cpt
Int_S=Int_S/cpt
print "Int Reg",Int_R
print "Int Sing",Int_S

#Reg
figure(1)
plt.plot(log(range(100,1000,10)),log(err1),"b")
plt.plot(log(range(100,1000,10)), -0.5*log(range(100,1000,10))+0.3, "r")
plt.title("Case 2D triangular regular $e^{-|x-y|^2}$|| pente = -0.5")
plt.xlabel("log(N)")
plt.ylabel("log(err)")
plt.show()

#Sing
figure(2)
plt.plot(log(range(100,1000,10)),log(err2),"b")
plt.plot(log(range(100,1000,10)), -0.5*log(range(100,1000,10))+0.25, "r")
plt.title("Case 2D triangular singular $1/|x-y|$|| pente = -0.5")
plt.xlabel("log(N)")
plt.ylabel("log(err)")
plt.show()

#time execution
figure(3)
plt.subplot(224)
plt.plot(range(100,N,10),temps,'-*r',linewidth=5)
plt.title("Triangle 2D")
plt.xlabel("N")
plt.ylabel("temps (sec)")
plt.show()