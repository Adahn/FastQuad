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



# returns points chosen randomly in the triangle A, B, C
#y=cst
def monteCarlo_Regul_Tria1(A,B,C,N):
    r1 = np.random.uniform(0, 1, N)
    r2 = np.random.uniform(0, 1, N)
    r1min=np.minimum(r1,r2)
    r2max=np.maximum(r1,r2)
    x=r1min*A[0]+(r2max-r1min)*B[0]+(1-r2max)*C[0]
    y=A[1]
    z=r1min*A[2]+(r2max-r1min)*B[2]+(1-r2max)*C[2]
    P=np.array([x,y,z])     
    return P

#z=cst
def monteCarlo_Regul_Tria2(A,B,C,N):
    r1 = np.random.uniform(0, 1, N)
    r2 = np.random.uniform(0, 1, N)
    r1min=np.minimum(r1,r2)
    r2max=np.maximum(r1,r2)
    x=r1min*A[0]+(r2max-r1min)*B[0]+(1-r2max)*C[0]
    y=r1min*A[1]+(r2max-r1min)*B[1]+(1-r2max)*C[1]
    z=A[2]
    P=np.array([x,y,z])     
    return P

#formule de Héron pour T1 y cst
def Aire_Tria_T1(A,B,C):
    a=np.sqrt((B[0]-C[0])**2+(B[2]-C[2])**2)
    b=np.sqrt((A[0]-C[0])**2+(A[2]-C[2])**2)
    c=np.sqrt((B[0]-A[0])**2+(B[2]-A[2])**2)
    P=a+b+c
    p=P/2
    S=np.sqrt(p*(p-a)*(p-b)*(p-c))
    return S

#formule de Héron pour T2 z cst
def Aire_Tria_T2(A,B,C):
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
    y1 = P1[1]
    z1 = P1[2]
    x2 = P2[0]
    y2 = P2[1]
    z2 = P2[2]
    #f_moy=sum(x1**2 + y1**2 + x2**2 + y2**2)/N
    #f_moy=sum(np.exp(-(np.sqrt((x1-x2)**2 + (y1-y2)**2))**2))
    f_moy=sum(np.exp(-((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)))/N
    f2=sum(np.exp(-((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2))**2)
    #Var=np.abs((1/N)*sum((x1**2 + y1**2 + x2**2 + y2**2)**2)-f_moy**2)
    Var=np.abs(f2-N*f_moy**2)/(N-1)
    err = np.sqrt(Var/N)*1.96/f_moy
    return err

def Int_reg(P1,P2,V):
    x1 = P1[0]
    y1 = P1[1]
    z1 = P1[2]
    
    x2 = P2[0]
    y2 = P2[1]
    z2 = P2[2]
    n = len(x1)
    f_moy=sum(np.exp(-((x1-x2)**2 + (y1-y2)**2 +(z1-z2)**2)))/n
    I=V*f_moy
    return I
    
# returns error for function f(x, y) = 1/|x-y|
def monteCarlo_Tria_sing(P1,P2,N):
    x1 = P1[0]
    y1 = P1[1]
    z1 = P1[2]
    x2 = P2[0]
    y2 = P2[1]
    z2 = P2[2]
    f_moy=sum(1/(np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)))/N
    f2=sum(1/(np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2))**2)
    Var=np.abs(f2-N*f_moy**2)/(N-1)
    err = np.sqrt(Var/N)*1.96/f_moy
    return err
    

def Int_sing(P1,P2,V):
    x1 = P1[0]
    y1 = P1[1]
    z1 = P1[2]
    
    x2 = P2[0]
    y2 = P2[1]
    z2 = P2[2]
    n = len(x1)
    
    f_moy=sum(1/(np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)))/n
    I=V*f_moy
    return I
# ----------------------------------------------------------------
# Calculate and plots
#-----------------------------------------------------------------

# number of points
N=1000
err1=np.array([])
err2 = np.array([])
for eps in np.array([0,1,1.99,2]):
    # triangle 1
    A1=np.array([2,0,1])
    B1=np.array([2.5,0,1])
    C1=np.array([2.2,0,5])
    print("Triangle( A1="+str(A1)+", B1="+str(B1)+", C1="+str(C1)+")")
    
    # triangle 2
    A2=np.array([2,eps,1])
    B2=np.array([2.5,eps,1])
    C2=np.array([2.2,eps,5])
    print("Triangle( A2="+str(A2)+", B2="+str(B2)+", C2="+str(C2)+")")
    
    # choose points randomly in triangle
    Point=monteCarlo_Regul_Tria1(A1,B1,C1,N)
    Point2=monteCarlo_Regul_Tria1(A2,B2,C2,N)
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
    err_Reg=np.array([])
    err_Sing = np.array([])
    
    for n in range(100,N,10):
        err1_=0
        err2_=0
        for k in range(10):
            P1=monteCarlo_Regul_Tria1(A1,B1,C1,n)
            P2=monteCarlo_Regul_Tria2(A2,B2,C2,n)
            err1_=err1_+ monteCarlo_Tria_regul(P1, P2, n)
            err2_=err2_+monteCarlo_Tria_sing(P1, P2, n)
        # calculate error
        err1_=err1_/10
        err2_=err2_/10
        err_Reg = np.append(err_Reg, err1_)
        err_Sing = np.append(err_Sing, err2_)
    err1=np.append(err_Reg,[err1])
    err2=np.append(err_Sing,[err2])
x=log(range(100,1000,10))

plt.figure(2)
plt.subplot(211)
plt.plot(x, log(err1[0:90]),label = str(2))
plt.plot(x, log(err1[90:180]),label=str(1))
plt.plot(x,log(err1[180:270]),label=str(0.01))
plt.plot(x,log(err1[270:360]),label=str(0))
plt.plot( x, -0.5*log(range(100,1000,10))+0.4,"r--")
plt.legend(bbox_to_anchor=(0.,1.02,1.,.102),loc=3,ncol=2,mode="expand",borderaxespad=0.)
plt.title("Pair of triangles regular polynom $exp^{-|x-y|^2}$|| pente = -0.5")
plt.xlabel("log(N)")
plt.ylabel("log(err)")
plt.subplot(212)
plt.plot(x, log(err2[0:90]),label = str(2))
plt.plot(x, log(err2[90:180]),label = str(1))
plt.plot(x,log(err2[180:270]),label = str(0.01))
plt.plot(x,log(err2[270:360]), label = str(0))
plt.plot(x, -0.5*log(range(100,1000,10))+0.2,"r--")
plt.legend(bbox_to_anchor=(0.,1.02,1.,.102),loc=3,ncol=2,mode="expand",borderaxespad=0.)
plt.title("Pair of triangles singular $1/(|x-y|)$ || pente = -0.5")
plt.xlabel("log(N)")
plt.ylabel("log(err)")

plt.show()
'''

# plot regular case
plt.figure()
plt.plot(log(range(100,N,10)),log(err_Reg),"b")
plt.plot(log(range(100,N,10)), -0.5*log(range(100,N,10))+0.5, "r")
plt.title("Pair of triangles regular case $exp^{-|x-y|^2}$|| pente = -0.5")
plt.xlabel("log(N)")
plt.ylabel("log(err)")

# plot singular case
plt.figure()
plt.plot(log(range(100,N,10)),log(err_Sing),"b")
plt.plot(log(range(100,N,10)), -0.5*log(range(100,N,10))+0.5, "r")
plt.title("Pair of triangles singular $1/(|x-y|)$ || pente = -0.5")
plt.xlabel("log(N)")
plt.ylabel("log(err)")
plt.show()
'''