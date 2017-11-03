# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 14:15:38 2016

@author: 3201955
"""

import random
import numpy as np
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
    
#random points in triangle ABC
def monteCarlo_triangle(A,B,C,N):
    r1 = np.random.uniform(0,1,N)
    r2 = np.random.uniform(0,1,N)
    x=(1-np.sqrt(r1))*A[0]+(np.sqrt(r1)*(1-r2))*B[0]+(r2*np.sqrt(r1))*C[0]
    y=(1-np.sqrt(r1))*A[1]+(np.sqrt(r1)*(1-r2))*B[1]+(r2*np.sqrt(r1))*C[1]
    P=np.array([x,y])     
    return P
    
def monteCarlo_quadrangle(A,B,C,D,N):
    P = monteCarlo_triangle(A,B,C,N)
    Q = monteCarlo_triangle(A,D,C,N)
    R = np.append(P,Q, axis=1)
    return R
    
A=np.array([1,0])
B=np.array([2,0])
C=np.array([2,0.9])
D=np.array([1,0.9])


def aire(A,B,C,D):
    V=(C[0]-D[0])*(C[1]-B[1])
    return V
    
#plt.figure()
#plt.scatter(A[0], A[1])
#plt.scatter(B[0], B[1])
#plt.scatter(C[0], C[1])
#plt.scatter(D[0],D[0])
#plt.scatter(Point[0],Point[1])
#plt.show()
    
def monteCarlo_Reg(P,N):
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
def monteCarlo_Sing(P,N):
    f_moy=sum(singular_function(P))/N
    f_moy2=sum(singular_function(P)**2)
    Var=np.abs(f_moy2-N*f_moy**2)/(N-1)
    err = np.sqrt(Var/N)*1.96/f_moy
   # print "f_err", err   
    return err

def Int_sing(Point1,V):
    x1, y1 = Point1
    N = len(x1)
    f_moy=sum(singular_function(Point1))/N
    I=V*f_moy
    return I

V=aire(A,B,C,D)
err1=np.array([])
err2=np.array([])
I_Reg=0
I_Sing=0
Int_R=0
Int_S=0
cpt=0
temps=np.array([])
for n in range(N1,N,pas):
    debut=time.time()
    cpt=cpt+1
    err1_=0
    err2_=0
    for k in range(10):
        Point=monteCarlo_quadrangle(A,B,C,D,n)
        err1_=err1_+monteCarlo_Reg(Point,n)
        err2_=err2_+monteCarlo_Sing(Point,n)
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
    tps=fin-debut
    temps=np.append(temps,tps)
    #print "temps", tps
Int_R=Int_R/cpt
Int_S=Int_S/cpt
print "Int Reg",Int_R
print "Int Sing",Int_S

#print size(err2),size(err1)
#Case reg
plt.plot(log(range(N1,N,pas)),log(err1),"b")
plt.plot(log(range(N1,N,pas)), -0.5*log(range(N1,N,pas))+0.07, "r")
plt.title("Case 2D quadrangle regular $exp(-|x-y|^{2})$")
plt.xlabel("log(N)")
plt.ylabel("log(err)")
plt.show()

#Case Sing
figure()
plt.plot(log(range(N1,N,pas)),log(err2),"b")
plt.plot(log(range(N1,N,pas)), -0.5*log(range(N1,N,pas))+0.1, "r")
plt.title("Case 2D quadrangle singular $1/|x-y|$")
plt.xlabel("log(N)")
plt.ylabel("log(err)")
plt.show()

figure(3)
plt.subplot(221)
plt.plot(range(N1,N,pas),temps,'-*r',linewidth=5)
plt.title("Rectangle 2D")
plt.xlabel("N")
plt.ylabel("temps (sec)")
plt.show()
