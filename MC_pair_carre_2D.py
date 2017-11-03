# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 14:15:38 2016

@author: 3201955
"""

import random
import numpy as np
import math
import matplotlib.pyplot as plt

N1=100 #borne inf
pas=10 #pas
N=1000

#regular function that we want to integrate 
def regular_function(Point1,Point2):
    x1, y1 = Point1
    x2, y2 = Point2
    f=np.exp(-(np.abs((x1-x2)**2+(y1-y2)**2))) # e^(-|x-y|²)
    return f

#singular function that we want to integrate
def singular_function(Point1,Point2):
    x1, y1 = Point1
    x2, y2 = Point2
    f=1/np.sqrt((x1-x2)**2+(y1-y2)**2) # 1/(|x-y|)
    return f
    
#random points in triangle ABC
def monteCarlo_triangle(A,B,C,N):
    r1 = np.random.uniform(0,1,N)
    r2 = np.random.uniform(0,1,N)
    x=(1-np.sqrt(r1))*A[0]+(np.sqrt(r1)*(1-r2))*B[0]+(r2*np.sqrt(r1))*C[0]
    y=(1-np.sqrt(r1))*A[1]+(np.sqrt(r1)*(1-r2))*B[1]+(r2*np.sqrt(r1))*C[1]
    #z=(1-np.sqrt(r1))*A[2]+(np.sqrt(r1)*(1-r2))*B[2]+(r2*np.sqrt(r1))*C[2]
    P=np.array([x,y])     
    return P
    
def monteCarlo_quadrangle(A,B,C,D,N):
    P = monteCarlo_triangle(A,B,C,N)
    Q = monteCarlo_triangle(A,D,C,N)
    R = np.append(P,Q, axis=1)
    return R

def aire(A,B,C,D):
    V=(C[0]-D[0])*(C[1]-B[1])
    return V
# funtion f(x,y)) = exp(-|x-y|^2)
def monteCarlo_Rec_regul(P1,P2,N):
    x1 = P1[0]
    N=len(x1)
    f_moy=sum(regular_function(Point1,Point2))/N
    f2=sum(regular_function(Point1,Point2)**2)
    Var=np.abs(f2-N*f_moy**2)/(N-1)
    err = np.sqrt(Var/N)*1.96/f_moy
    return err

def Int_reg(Point1,Point2,V):
    x1, y1 = Point1
    n = len(x1)
    f_moy=sum(regular_function(Point1,Point2))/n
    I=V*f_moy
    return I
    
# returns error for function f(x, y) = 1/|x-y|
def monteCarlo_Rec_sing(P1,P2,N):
    x1 = P1[0]
    N=len(x1)
    f_moy=sum(singular_function(Point1,Point2))/N
    f2=sum(singular_function(Point1,Point2)**2)
    Var=np.abs(f2-N*f_moy**2)/(N-1)
    err = np.sqrt(Var/N)*1.96/f_moy
    return err

def Int_sing(Point1,Point2,V):
    x1, y1 = Point1
    n = len(x1)
    
    f_moy=sum(singular_function(Point1,Point2))/n
    I=V*f_moy
    return I
############################################################
############################################################

# rectangle 1
A1=np.array([1,1])
B1=np.array([4,1])
C1=np.array([4,3])
D1=np.array([1,3])
print("Rectangle1:  ( A1="+str(A1)+", B1="+str(B1)+", C1="+str(C1)+", D1="+str(D1)+")")

# rectangle 2
A2=np.array([1,-2])
D2=np.array([1,2])
C2=np.array([4,2])
B2=np.array([4,-2])
print("Rectangle2:  ( A2="+str(A2)+", B2="+str(B2)+", C2="+str(C2)+", D2="+str(D2)+")")

#
#plt.figure()
#plt.scatter(A1[0], A1[1])
#plt.scatter(B1[0], B1[1])
#plt.scatter(C1[0], C1[1])
#plt.scatter(D1[0], D1[1])
#plt.scatter(Point1[0],Point1[1])
#plt.scatter(A2[0], A2[1],c='r')
#plt.scatter(B2[0], B2[1],c='r')
#plt.scatter(C2[0], C2[1],c='r')
#plt.scatter(D2[0], D2[1],c='r')
#plt.scatter(Point2[0],Point2[1],c='r')

############################################################

err_Reg=np.array([])
err_Sing = np.array([])

V=aire(A1,B1,C1,D1)*aire(A2,B2,C2,D2)
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
        Point1=monteCarlo_quadrangle(A1,B1,C1,D1,n)
        #print(size(Point1))
        Point2=monteCarlo_quadrangle(A2,B2,C2,D2,n)
        #print(size(Point2))
        err1_=err1_+monteCarlo_Rec_regul(Point1,Point2,n)
        err2_=err2_+monteCarlo_Rec_sing(Point1,Point2,n)
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
plt.plot(log(range(N1,N,pas)), -0.5*log(range(N1,N,pas))+0.4, "r")
plt.title("Pair of rectangles regular polynom $exp^{-|x-y|^2}$|| pente = -0.5")
plt.xlabel("log(N)")
plt.ylabel("log(err)")

# plot singular case
plt.figure()
plt.plot(log(range(N1,N,pas)),log(err2),"b")
plt.plot(log(range(N1,N,pas)), -0.5*log(range(N1,N,pas))+0.2, "r")
plt.title("Pair of rectangles singular $1/(|x-y|)$ || pente = -0.5")
plt.xlabel("log(N)")
plt.ylabel("log(err)")


