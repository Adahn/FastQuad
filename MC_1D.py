# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 15:06:20 2017

@author: 3201955
"""

import random
import numpy as np
import math
import matplotlib.pyplot as plt
import time
N = 1000
a = -1
b = 1
'''
# Case 1 : regular function constant f(x) = 4
def monteCarlo_Regul_Const(a, b, N):
    randValues = np.random.uniform(a, b, N)
    f_moy = (4*N)/N
    f_int = (b-a)*f_moy
    f_moy2 = (4**2*N)/N
    
    error = np.sqrt((f_moy2-f_moy**2)/N) / f_int
    return error
   '''
# Case 1 : regular function polynom f(x) = x^2
def monteCarlo_Regul_poly(a, b, N):
    X = np.random.uniform(a, b, N)
    f_moy = sum(X**2)/N
    f_int = (b-a)*f_moy
    f_moy2 = sum((X**2)**2)
    error = (b-a)*np.sqrt(((f_moy2 - N*f_moy**2)/(N-1))/N)*1.96 / f_int
    #print "int_reg",f_int
    return error
    
#case 2: regular function exp(x)
def monteCarlo_Regul_exp(a, b, N):
    X = np.random.uniform(a, b, N)
    f_moy = sum(np.exp(X))/N
    f_int = (b-a)*f_moy
    f_moy2 = sum((np.exp(X))**2)
    error = (b-a)*np.sqrt(((f_moy2 - N*f_moy**2)/(N-1))/N)*1.96 / f_int
    #print "int_reg",f_int
    return error
def Int_reg(Point1,V):
    x1 = Point1
    n = len(x1)
    f_moy = sum(np.exp(x1))/n
    I=V*f_moy
    return I

# Case 3 : singular function f(x)  =|x|
def monteCarlo_Sing(a, b, N):
    randValues = np.random.uniform(a, b, N)
    f_moy = sum(np.abs(randValues))/N
    f_int = (b-a)*f_moy
    f_moy2 = sum(np.abs(randValues)**2)
    Var=(f_moy2-N*f_moy**2)/(N-1)
    error = np.sqrt(Var/N)*1.96/f_moy
   # print "int_sing",f_int
    return error

def Int_sing(Point1,V):
    x1= Point1
    n = len(x1)
    f_moy = sum(np.abs(x1))/n
    I=V*f_moy
    return I
#--------------------------------------------------
V=b-a
temps=np.array([])
errVec1 = np.array([])
errVec2 = np.array([])
I_Reg=0
I_Sing=0
Int_R=0
Int_S=0
Point=np.array([])
for n in range(100,1000,10):
    debut=time.time()
    err1_=0
    err2_=0
    for k in range(10):
        Point=np.random.uniform(a, b, n)
        err1_ = err1_+monteCarlo_Regul_exp(a, b, n)
        err2_ = err2_+monteCarlo_Sing(a, b, n)
        I_Reg=I_Reg+Int_reg(Point,V)
        I_Sing=I_Sing+Int_sing(Point,V)
    I_Reg=I_Reg/10
    I_Sing=I_Sing/10
    Int_R=Int_R+I_Reg
    Int_S=Int_S+I_Sing
    err1_=err1_/10
    err2_=err2_/10
    errVec1 = np.append(errVec1, err1_)
    errVec2 = np.append(errVec2, err2_)
    fin=time.time()
    tps=fin-debut
    temps=np.append(temps,tps)
Int_R=Int_R/cpt
Int_S=Int_S/cpt
print "Int Reg",Int_R
print "Int Sing",Int_S



#time exec
figure(3)
subplot(222)
plt.plot(range(100,1000,10),temps,"-*r",linewidth=5)
plt.title("1D")
plt.xlabel("N")
plt.ylabel("temps (sec)")
plt.show()

#figure(1)
#errVec2 = np.array([])
#print("Case regular case")
#for n in range(100,1000,10):
#    err_=0
#    for k in range(10):
#        err_ = err_+monteCarlo_Regul_poly(a, b, n)
#    err_=err_/10
#    errVec2 = np.append(errVec2, err_)

figure(1)
plt.plot(np.log(range(100,1000,10)), np.log(errVec1))
plt.plot(np.log(range(100,1000,10)), -0.5*np.log(range(100,1000,10))+0.55, "r")
plt.title("Case 1D regular case $exp^{x}$ || pente = -0.5")
plt.xlabel("log(N)")
plt.ylabel("log(err)")
plt.show()

figure(2)
plt.plot(np.log(range(100,1000,10)), np.log(errVec2))
plt.plot(np.log(range(100,1000,10)), -0.5*np.log(range(100,1000,10))+0.55, "r")
plt.title("Case 1D singular case $|x|$ || pente = -0.5")
plt.xlabel("log(N)")
plt.ylabel("log(err)")
plt.show()