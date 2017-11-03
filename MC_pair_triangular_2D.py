# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 14:09:50 2016

@author: Adrian Ahne

Calculate integral of function over pair of triangular
"""


import random 
import numpy as np
import math
import matplotlib.pyplot as plt
from numpy import linalg as LA
from numpy.linalg import norm

N1=100 #borne inf
pas=10 #pas
N=1000
#regular function that we want to integrate 
def regular_function(Point1,Point2):
    x1,y1=Point1
    x2,y2=Point2
    f=np.exp(-(np.abs((x1-x2)**2+(y1-y2)**2))) # e^(-|x-y|²)
    return f

#singular function that we want to integrate
def singular_function(Point1,Point2):
    x1,y1=Point1
    x2,y2=Point2
    f=1/sqrt(np.abs((x1-x2)**2+(y1-y2)**2)) # 1/(|x-y|)
    return f
    

# returns points chosen randomly in the triangle A, B, C
#y=cst
def monteCarlo_Tria(A,B,C,N):
    r1 = np.random.uniform(0, 1, N)
    r2 = np.random.uniform(0, 1, N)
    r1min=np.minimum(r1,r2)
    r2max=np.maximum(r1,r2)
    x=r1min*A[0]+(r2max-r1min)*B[0]+(1-r2max)*C[0]
    y=r1min*A[1]+(r2max-r1min)*B[1]+(1-r2max)*C[1]
    P=np.array([x,y])     
    return P

#formule de Héron pour T2 z cst
def Aire_Tria(A,B,C):
    a=np.sqrt((B[0]-C[0])**2+(B[1]-C[1])**2)
    b=np.sqrt((A[0]-C[0])**2+(A[1]-C[1])**2)
    c=np.sqrt((B[0]-A[0])**2+(B[1]-A[1])**2)
    P=a+b+c
    p=P/2
    S=np.sqrt(p*(p-a)*(p-b)*(p-c))
    return S
    
# funtion f(x,y)) = exp(-|x-y|^2)
def monteCarlo_Tria_regul(P1,P2,N):
    x1 = P1[0]

    N=len(x1)
    f_moy=sum(regular_function(P1,P2))/N
    f2=sum(regular_function(P1,P2)**2)
    Var=np.abs(f2-N*f_moy**2)/(N-1)
    err = np.sqrt(Var/N)*1.96/f_moy
    return err

def Int_reg(P1,P2,V):
    x1, y1 = P1
    N = len(x1)
    
    f_moy=sum(regular_function(P1,P2))/N
    I=V*f_moy
    return I
# returns error for function f(x, y) = 1/|x-y|
def monteCarlo_Tria_sing(P1,P2,N):
    x1 = P1[0]
    N=len(x1)
    f_moy=sum(singular_function(P1,P2))/N
    f2=sum(singular_function(P1,P2)**2)
    Var=np.abs(f2-N*f_moy**2)/(N-1)
    err = np.sqrt(Var/N)*1.96/f_moy
    return err
    
def Int_sing(P1,P2,V):
    x1, y1 = Point1
    N = len(x1)
    
    f_moy = sum(singular_function(P1,P2))/N
    I=V*f_moy
    return I
    
# ----------------------------------------------------------------
# Calculate and plots
#-----------------------------------------------------------------

# triangle 1
A1=np.array([2,1])
B1=np.array([2.5,1])
C1=np.array([2.2,5])
print("Triangle( A1="+str(A1)+", B1="+str(B1)+", C1="+str(C1)+")")

# triangle 2
A2=np.array([2,-1])
B2=np.array([2.2,1])
C2=np.array([2.5,-1])
print("Triangle( A2="+str(A2)+", B2="+str(B2)+", C2="+str(C2)+")")

# choose points randomly in triangle
#Point=monteCarlo_Regul_Tria(A1,B1,C1,N)
#Point2=monteCarlo_Regul_Tria(A2,B2,C2,N)
#
#plt.figure()
#plt.scatter(A1[0], A1[2])
#plt.scatter(B1[0], B1[2])
#plt.scatter(C1[0], C1[2])
#plt.scatter(Point[0],Point[2])
#plt.scatter(A2[0], A2[1],c='r')
#plt.scatter(B2[0], B2[1],c='r')
#plt.scatter(C2[0], C2[1],c='r')
#plt.scatter(Point2[0],Point2[1],c='r')
#plt.show()
#print "Aire T1:",Aire_Tria_T1(A1,B1,C1)
#print "Aire T2:",Aire_Tria_T2(A2,B2,C2)
#Aire=Aire_Tria_T1(A1,B1,C1)*Aire_Tria_T2(A2,B2,C2)
#print "Aire:", Aire

V=Aire_Tria(A1,B1,C1)*Aire_Tria(A2,B2,C2)

err1=np.array([])
err2=np.array([])

I_Reg=0
I_Sing=0
Int_R=0
Int_S=0
cpt=0
for n in range(N1,N,pas):
    cpt=cpt+1
    err1_=0
    err2_=0
    for k in range(10):
        Point1=monteCarlo_Tria(A1,B1,C1,n)
        #print(size(Point1))
        Point2=monteCarlo_Tria(A2,B2,C2,n)
        #print(size(Point2))
        err1_=err1_+monteCarlo_Tria_regul(Point1,Point2,n)
        err2_=err2_+monteCarlo_Tria_sing(Point1,Point2,n)
        I_Reg=I_Reg+Int_reg(Point1,Point2,V)
        I_Sing=I_Sing+Int_sing(Point1,Point2,V)
    I_Reg=I_Reg/10
    I_Sing=I_Sing/10
    Int_R=Int_R+I_Reg
    Int_S=Int_S+I_Sing
    err1_=err1_/10
    err2_=err2_/10
    err1=np.append(err1,err1_)
    err2=np.append(err2,err2_)
Int_R=Int_R/cpt
Int_S=Int_S/cpt
print "Int Reg",Int_R
print "Int Sing",Int_S


# plot regular case
plt.figure()
plt.plot(log(range(N1,N,pas)),log(err1),"b")
plt.plot(log(range(N1,N,pas)), -0.5*log(range(N1,N,pas))+1.4, "r")
plt.title("Pair of triangles regular case $exp^{-|x-y|^2}$|| pente = -0.5")
plt.xlabel("log(N)")
plt.ylabel("log(err)")

# plot singular case
plt.figure()
plt.plot(log(range(N1,N,pas)),log(err2),"b")
plt.plot(log(range(N1,N,pas)), -0.5*log(range(N1,N,pas))+0.2, "r")
plt.title("Pair of triangles singular $1/(|x-y|)$ || pente = -0.5")
plt.xlabel("log(N)")
plt.ylabel("log(err)")
