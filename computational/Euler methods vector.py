# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 11:31:09 2024

@author: Usr
"""
#Euler differential equation methods using vectors
import scipy as sci
import numpy as np
import matplotlib.pyplot as plt

alpha = 0.1
omega = 2
t = np.linspace(0,50, 10000)
x0 = 1
v0 = 0

delta_t = (t[1]-t[0])
print(delta_t)

def Euler_forward(x0,v0,t,alpha,omega,delta_t):
    x = np.zeros(len(t))
    v = np.zeros(len(t))
    x[0],v[0] = x0,v0
    for i in range(0,len(t)-1):
        x[i+1] = x[i] +delta_t*v[i]
        v[i+1] = v[i] + delta_t*(-omega**2*x[i] -alpha*v[i])
    return x,v

def Euler_backward(x0,v0,t,alpha,omega,delta_t):
    x = np.zeros(len(t))
    v = np.zeros(len(t))
    x[0],v[0] = x0,v0
    for j in range(0,len(t)-1):
        v[j+1] = (v[j] - (omega**2*x[j]*delta_t))/(1+alpha*delta_t)
        x[j+1] = (x[j] + delta_t*(v[j+1]))
    return x,v

def Trapezoidal(x0,v0,t,alpha,omega,delta_t):
    v = np.zeros(len(t))
    x = np.zeros(len(t))
    v[0], x[0] = v0, x0
    delta_t = (t[-1]-t[0])/len(t)
    #hay que utilizar euler forward para predecir los valores v[i+1] y x[i+1]
    for m in range(0,len(t)-1):
        v_antes = v[m] + delta_t*(-omega**2*x[m] -alpha*v[m])
        x[m+1] = x[m]+ (delta_t/2*(v[m]+v_antes))
        v[m+1] = v[m] + (delta_t/2)*(-omega**2*x[m+1]-alpha*v_antes-omega**2*x[m]-alpha*v[m])
    return x,v

def Exacta(x0,v0,t,alpha,omega,delta_t):
    Omega = np.sqrt(omega**2 - alpha**2 / 4)
    C1 = x0
    C2 = (v0 + alpha * x0 / 2) / Omega
    return np.exp(-alpha * t / 2) * (C1 * np.cos(Omega * t) + C2 * np.sin(Omega * t))

x1,v1= Euler_forward(x0,v0,t,alpha, omega,delta_t)
x2,v2 = Euler_backward(x0,v0,t,alpha, omega,delta_t)
x3,v3 = Trapezoidal(x0,v0, t,alpha, omega,delta_t)
plt.plot(t,x1)
plt.plot(t,x2)
plt.plot(t,x3)
plt.show()
    

  