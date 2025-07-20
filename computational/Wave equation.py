# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 12:50:55 2024

@author: Usr
"""
#Wave equation solution 1D
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, Animation
from matplotlib.animation import PillowWriter
from scipy import sparse

b = 0.2 #2k/gamma
a = -2
c = 0.1
x = np.linspace(0,2,500)
dx = 0.01
t = np.linspace(0,40,500)
dt = 0.01
#subamortiguado a = -2, dt = 0.001, dx = 0.1
#sobreamortiguado a = -1, dt y dx igual

n = [1,2,3,4,5]

CFL = c * dt / dx
if CFL > 1:
    print(f"Criterio CFL no satisfecho. CFL = {CFL}")

def EcuacionOnda1D(t,x,dt,dx,a,b,c,n):
    
    alpha = (c**2 * dt**2)/dx**2
    beta = 2+ (a*dt**2)
    d = (-2*(alpha)+beta)*np.ones(len(x))

    u = alpha*np.ones(len(x)-1)
    A = sparse.diags([u,d,u],[1,0,-1],  shape = (len(x),len(x))).toarray()   

    gamma= -1 +(b*dt/2)
    B = sparse.eye(len(x))*gamma
    u = np.sin(n*np.pi*x)
    soluciones = [np.copy(u)]
    u_antes = np.copy(u)
    for j in range(1,len(x)):
        u_3 = ((A @ u) +B@u_antes)/(1+(b*dt/2))
        u_3[0] = u_3[-1] = 0
        
        soluciones.append(np.copy(u_3))
        u_antes = np.copy(u)
        u = np.copy(u_3)
    return np.array(soluciones)

for i in n:
    soluciones = EcuacionOnda1D(t,x,dt,dx,a,b,c,i)


    fig, ax = plt.subplots()
    line, = ax.plot(x, soluciones[0])
    ax.set_title('Ecuación onda 1D Diferencias centradas')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('u(x,t)')
    def animate(i):
        
        line.set_ydata(soluciones[i])  
    
    
    
        return line,

    ani = FuncAnimation(
        fig, animate, frames = len(x), interval=1, blit=True)
    ani.save(f'Diferencias centradas{i}.gif', writer=PillowWriter(fps = 50))
    plt.show()
plt.close()


n = [1,2,3,4,5]
def EcuacionOndaCN(t,x,dt,dx,a,b,c,n):
    alpha = (c**2 * dt**2)/dx**2
    beta = 2 + (a*dt**2) +(b*dt/2) 
    d = (-2*(alpha)+beta)*np.ones(len(x))

    u = alpha*np.ones(len(x)-1)
    A = sparse.diags([u,d,u],[1,0,-1],  shape = (len(x),len(x))).toarray()   
    
    
    B = sparse.eye(len(x))
    u = np.sin(n*np.pi*x)
    solucionesCN = [np.copy(u)]
    u_antes = np.copy(u)

    
    for j in range(1,len(x)):
        u_3 = ((A @ u) - B@u_antes)/(1+(b*dt/2))
        u_3[0] = u_3[-1] = 0
        
        solucionesCN.append(np.copy(u_3))
        u_antes = np.copy(u)
        u = np.copy(u_3)
        
   
    return np.array(solucionesCN)

for j in n:
    
    solucionesCN = EcuacionOndaCN(t,x,dt,dx,a,b,c,j)



    fig, ax = plt.subplots()
    line, = ax.plot(x, solucionesCN[0])
    ax.set_title('Ecuación onda 1D Crank-Nicolson')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('u(x,t)')

    def animate(j):
        
        line.set_ydata(solucionesCN[j])  
    
        return line,

    ani = FuncAnimation(
        fig, animate,frames = len(x), interval=10, blit=True)
    ani.save(f'Crank-Nicolson {j}.gif', writer=PillowWriter(fps = 50))
    plt.show()
plt.close()