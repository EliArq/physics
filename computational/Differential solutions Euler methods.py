# -*- coding: utf-8 -*-

#Differential equation solutions using Euler methods using matrix solve
import numpy as np
import matplotlib.pyplot as plt

omega = 2
alpha = 0.1

A = np.array([[0,1],[-omega**2,-alpha]])

t = np.linspace(0,50, 20000)
x0 = 1
v0 = 0

dt = (t[1]-t[0])


def Euler_forward(v0,x0,A,t,dt):
    u = np.zeros((2, len(t))) 
    u[0, 0] = x0  
    u[1, 0] = v0  
    I = np.eye(2)
    A1 = I +dt*A
    for i in range(0, len(t) - 1):
        u[:, i+1] = (I + dt*A)@u[:,i]
    return u,A1

def Euler_backward(v0,x0,A,t,dt):
    u = np.zeros((2, len(t))) 
    u[0, 0] = x0  
    u[1, 0] = v0  
    I = np.eye(2)
    B = (I - A*dt)
    for i in range(0, len(t) - 1):
       u[:,i+1] = np.linalg.solve((I - A*dt),u[:,i])
    return u,B

def Trapezoidal(v0,x0,A,t,dt):
    u = np.zeros((2, len(t))) 
    u[0, 0] = x0  
    u[1, 0] = v0  
    I = np.eye(2)
    C = (I +dt/2*A)@(I -dt/2*A)
    for i in range(0, len(t) - 1):
        derecho = (I+ dt/2*A) @ u[:,i]
        izquierdo = (I - dt/2*A)
        u[:,i+1] = np.linalg.solve(izquierdo,derecho)
    return u,C

def Exacta(v0,x0,A,t,dt):
    Omega = np.sqrt(omega**2-(alpha/2)**2)
    return np.exp(-alpha*t/2)*(x0*np.cos(Omega*t)+(v0+alpha*x0/2)*np.sin(Omega*t)/Omega)

u1,A1 = Euler_forward(v0,x0,A,t,dt)
u2,B = Euler_backward(v0,x0,A,t,dt)
u3,C = Trapezoidal(v0,x0,A,t,dt)
exacta = Exacta(v0,x0,A,t,dt)

#error of methods
errorEf = np.abs(u1[0,:]-exacta)
errorEb = np.abs(u2[0,:]-exacta)
errorT = np.abs(u3[0,:]-exacta)

plt.plot(t, u1[0, :], label='Euler Forward')  
plt.plot(t, u2[0, :], label='Euler Backward') 
plt.plot(t, u3[0, :], label='Crank-Nicholson') 
plt.plot(t,exacta,label = 'Solución exacta', linestyle= '--')
plt.xlabel('Tiempo t(s)')
plt.ylabel('Posición x(m)')
plt.legend()
plt.show()


plt.plot(t,errorEf, label = 'Euler Forward error', linestyle= '--')  
plt.plot(t, errorEb, label='Euler Backward error', linestyle= '--') 
plt.plot(t, errorT, label='Crank-Nicholson error', linestyle= '--') 
plt.xlabel('Paso de tiempo t(s)')
plt.ylabel('Incremento del error (m)')
plt.legend()
plt.show()

#stability
tiempoerror = np.array([1e-4, 1e-3, 5e-3, 1e-2, 2e-2])
erroresEf = np.zeros(len(tiempoerror))
erroresEb = np.zeros(len(tiempoerror))
erroresT = np.zeros(len(tiempoerror)) 
for i in range(len(tiempoerror)):
    u1,A1 = Euler_forward(v0,x0,A,t,tiempoerror[i])
    u2,B = Euler_backward(v0,x0,A,t,tiempoerror[i])
    u3,C = Trapezoidal(v0,x0,A,t,tiempoerror[i])
    exacta = Exacta(v0,x0,A,t,tiempoerror[i])
    errorEf = np.abs(u1[0,:]-exacta)
    errorEb = np.abs(u2[0,:]-exacta)
    errorT = np.abs(u3[0,:]-exacta)
    erroresEf[i] = max(np.abs(np.linalg.eigvals(A1)))
    erroresEb[i] = max(np.abs(np.linalg.eigvals(B)))
    erroresT[i] = max(np.abs(np.linalg.eigvals(C)))
    
    
plt.plot(tiempoerror,erroresEf, label = 'Estabilidad Euler forward')
plt.plot(tiempoerror,erroresEb, label = 'Estabilidad Euler backward')
plt.plot(tiempoerror,erroresT, label = 'Estabilidad Crank-Nicholson')
plt.xscale('log')
plt.xlabel('Paso de tiempo t(s)')
plt.ylabel('Valor de $| \lambda | (s)^(-1) $ ')
plt.legend()
plt.show()
print('Como se puede observar, tanto Euler backward como Crank-Nicholson son más estables con pasos de tiempo mayores, mientras que Euler forward diverge rápidamente con pasos mayores de tiempo')