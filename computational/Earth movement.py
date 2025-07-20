# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#Equation solves for Earth dynamic movement
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sci
#datos
G = 6.6738e-11
M = 1.9891e30
m = 5.9722e24

h = 3600
#h = 3600*60 
#a partir de este paso de tiempo runge-kutta no conserva la energía, es h considerablemente grande para la estabilidad del método
nt = 5*365*24*3600
dt = 3600
t = np.arange(0,nt,dt)
N = len(t)
Ec = []
Ep = []
def f(r):
    
    dist = np.sqrt((r[0]**2) + (r[1]**2))
   
    ax = -G*M*r[0]/dist**3
    ay = -G*M*r[1]/dist**3
   
    return np.array([ax,ay])

r = np.zeros((N,2))

r[0] = np.array([1.4719e11,0])
Ep.append(-G*M*m/np.sqrt(r[0][0]**2+r[0][1]**2))
v = np.zeros((N,2))

v[0] = np.array([0,3.0287e4])
Ec.append(0.5*m*(v[0][0]**2+v[0][1]**2))
v[1] = v[0] + h*f(r[0])*0.5
Ec.append(0.5*m*(v[1][0]**2+v[1][1]**2))
r[1] = r[0] + (h*v[1])
Ep.append(-G*M*m/np.sqrt(r[1][0]**2+r[1][1]**2))
for i in range(1,N-1):
    
    r[i+1] = r[i] + (h*v[i])
    
    k = h*f(r[i+1])
    v[i+1] = v[i] + h*f(r[i+1])
    Ec.append(0.5*m*(v[i+1][0]**2+v[i+1][1]**2))
    Ep.append(-G*M*m/np.sqrt(r[i+1][0]**2+r[i+1][1]**2))

plt.title('Órbita de la Tierra alrededor del Sol')
plt.plot(r[:,0],r[:,1])
plt.ylabel('r (m)')
plt.xlabel('t (s)')

plt.show()

Ec = np.array(Ec)
Ep = np.array(Ep)
plt.title('Energía del sistema')
plt.plot(t,Ec, label = 'Energía cinética')
plt.plot(t,Ep, label = 'Energía potencial')
plt.plot(t,Ec+Ep, label = 'Energía total')
plt.xlabel('t (s)')
plt.ylabel('E (J)')
plt.legend(loc = 'best')
plt.show()
print('Como se puede observar, la energía es constante y la órbita elíptica, \n aunque si se observa aumentando el zoom se puede ver que la órbita no pasa por el mismo punto siempre \n'
      ' y con la energía pasa algo similar, no es exactamente constante en todo momento.')


#runge-kutta 4

t = np.arange(0,nt,dt)
N = len(t)
r = np.zeros((N,2))
Ec = []
Ep = []
r[0] = np.array([1.4719e11,0])
Ep.append(-G*M*m/np.sqrt(r[0][0]**2+r[0][1]**2))
v = np.zeros((N,2))

v[0] = np.array([0,3.0287e4])
Ec.append(0.5*m*(v[0][0]**2+v[0][1]**2))

for i in range(0,N-1):
    k1 = v[i]
    k2 = v[i] +0.5*h*f(r[i])
    k3 = v[i] + 0.5*h*f(r[i] + 0.5*h*k2)
    k4 = v[i] + h*f(r[i] + h*k3)
    q1 = f(r[i])
    q2 = f(r[i] + 0.5*h*k1)
    q3 = f(r[i] + +0.5*h*k2)
    q4 = f(r[i] +h*k3)
    v[i+1] = v[i] + (h/6)*(q1+2*q2+2*q3+q4)
    r[i+1] = r[i] + (h/6)*(k1+2*k2+2*k3+k4)
    Ec.append(0.5*m*(v[i+1][0]**2+v[i+1][1]**2))
    Ep.append(-G*M*m/np.sqrt(r[i+1][0]**2+r[i+1][1]**2))

plt.title('Órbita de la Tierra alrededor del Sol Runge-Kutta 4')
plt.plot(r[:,0],r[:,1])
plt.ylabel('r (m)')
plt.xlabel('t (s)')

plt.show()

Ec = np.array(Ec)
Ep = np.array(Ep)
plt.title('Energía del sistema Runge-Kutta 4')
plt.plot(t,Ec, label = 'Energía cinética')
plt.plot(t,Ep, label = 'Energía potencial')
plt.plot(t,Ec+Ep, label = 'Energía total')
plt.xlabel('t (s)')
plt.ylabel('E (J)')
plt.legend(loc = 'best')
plt.show()


#runge-kutta 3
t = np.arange(0,nt,dt)
N = len(t)
r = np.zeros((N,2))
Ec = []
Ep = []
r[0] = np.array([1.4719e11,0])
Ep.append(-G*M*m/np.sqrt(r[0][0]**2+r[0][1]**2))
v = np.zeros((N,2))

v[0] = np.array([0,3.0287e4])
Ec.append(0.5*m*(v[0][0]**2+v[0][1]**2))

for i in range(0,N-1):
    k1 = v[i]
    k2 = v[i] +0.5*h*f(r[i] + 0.5*h*k1)
    k3 = v[i] + h*f(r[i] + h*k2)
    q1 = f(r[i])
    q2 = f(r[i] + 0.5*h*k1)
    q3 = f(r[i] + h*k2)
    
    v[i+1] = v[i] + (h/6)*(q1+4*q2+q3)
    r[i+1] = r[i] + (h/6)*(k1+4*k2+k3)
    Ec.append(0.5*m*(v[i+1][0]**2+v[i+1][1]**2))
    Ep.append(-G*M*m/np.sqrt(r[i+1][0]**2+r[i+1][1]**2))

plt.title('Órbita de la Tierra alrededor del Sol Runge-Kutta 3')
plt.plot(r[:,0],r[:,1])
plt.ylabel('r (m)')
plt.xlabel('t (s)')

plt.show()

Ec = np.array(Ec)
Ep = np.array(Ep)
plt.title('Energía del sistema Runge-Kutta 3')
plt.plot(t,Ec, label = 'Energía cinética')
plt.plot(t,Ep, label = 'Energía potencial')
plt.plot(t,Ec+Ep, label = 'Energía total')
plt.xlabel('t (s)')
plt.ylabel('E (J)')
plt.legend(loc = 'best')
plt.show()
print('Con Runge-Kutta 4 y 3 tenemos el mismo caso a Verlet, tanto para la órbita como la energía.')

#scipy solve_ivp
t_sim = (0,nt)
r = np.zeros((N,2))
t = np.arange(0,nt,dt)
Ec = []
Ep = []
r[0] = np.array([1.4719e11,0])
Ep.append(-G*M*m/np.sqrt(r[0][0]**2+r[0][1]**2))
v = np.zeros((N,2))
v[0] = np.array([0,3.0287e4])
Ec.append(0.5*m*(v[0][0]**2+v[0][1]**2))

def f(t,y):
    r = [y[0],y[1]]
    v = [y[2],y[3]]
    dist = np.sqrt(r[0]**2+r[1]**2)
    dvdt_x = -G*M*r[0]/dist**3
    dvdt_y = -G*M*r[1]/dist**3
    drdt_x = v[0]
    drdt_y = v[1]
    
    
    return [drdt_x,drdt_y,dvdt_x,dvdt_y]


y0 = np.concatenate([r[0],v[0]])
sol = sci.solve_ivp(f,t_sim,y0,t_eval = t,method = 'RK45')
rx = sol.y[0]
ry = sol.y[1]
plt.title('Órbita de la Tierra alrededor del Sol ivp')
plt.ylabel('r (m)')
plt.xlabel('t (s)')
plt.plot(rx,ry)
plt.show()
print('Como se observa, la librería que implementa scipy no satisface una órbita elíptica cerrada, \n'
      ' mostrando una órbita en espiral que no corresponde al movimiento observado de la Tierra alrededor del Sol.')

vx = sol.y[2]
vy = sol.y[3]

#Energy plot
for j in range(1,len(t)):
 
    Ec.append(0.5*m*(vx[j]**2+vy[j]**2))
    Ep.append(-G*M*m/np.sqrt(rx[j]**2+ry[j]**2))

Ec = np.array(Ec)
Ep = np.array(Ep)
plt.title('Energía sistema ivp')
plt.plot(t,Ec, label= 'Energía cinética')
plt.plot(t,Ep, label= 'Energía potencial')
plt.plot(t, Ec+Ep, label = 'Energía total')
plt.xlabel('t (s)')
plt.ylabel('E (J)')
plt.legend(loc = 'best')
plt.show()
print('Respecto a la energía, también se observa una diferencia respecto a los otros métodos, pues, como \n'
      'ya se representaba en la órbita, la energía no se mantiene constante, de hecho, decrece según pasa el tiempo,\n'
      ' la Tierra caería en el campo gravitatorio del Sol y sería absorbida por él.')