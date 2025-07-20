# -*- coding: utf-8 -*-

#Diffusion equation 1D using matrix
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from matplotlib.animation import FuncAnimation


x = np.linspace(-1,1,100)
y = np.linspace(-1,1,100)
D = 0.01
dx = (x[1]-x[0])
dt = 0.001
T = 100*np.exp(-20*x**2)
T_desp = np.copy(T)
estabilidad = dx**2/(2*D)

def Ecuaciondifusion(T,D,x,dx,dt):
    d = 2*np.ones(len(x))
    u = np.ones(len(x)-1)
    A = sparse.diags([d,u,u], [0,-1,1], shape = (len(x),len(x))).toarray()
    alpha = D*dt/dx**2
    I =sparse.eye(len(x)).toarray()
    M = I + alpha*A
    M[0,:] = 0
    M[-1,:] = 0
    M[:,0] = 0
    M[:,-1] = 0
    derivadas = M@T
    return derivadas

sol = Ecuaciondifusion(T,D,x,dx,dt)
fig, ax = plt.subplots()
line, = ax.plot(x, sol)
ax.set_title('Ecuación difusión 1D FTCS')
ax.set_xlabel('x (m)')
ax.set_ylabel('T')
def animate(i):
    global T_desp
    T_desp = Ecuaciondifusion(T_desp,D,x,dx,dt)
    line.set_ydata(T_desp)  
    return line,


ani = FuncAnimation(
    fig, animate, interval=60, blit=False)

plt.show()


#b

D = 0.1
dx = (x[1]-x[0])
T = 100*np.exp(-20*x**2)
T_desp = np.copy(T)
def Trapezoidal(T,D,x,dx,dt):
    alpha = D*dt/dx**2
    d = 1-2*alpha*np.ones(len(x))
    u = alpha*np.ones(len(x)-1)
    A = sparse.diags([d,u,u], [0,-1,1], shape = (len(x),len(x))).toarray()
    o = 1+2*alpha*np.ones(len(x))
    B = sparse.diags([o,-u,-u],[0,1,-1], shape = (len(x),len(x))).toarray()
    B[0] = np.zeros(len(x))
    B[0,0] = 1
    B[-1,-2] = -1
    B[-1,-1] = 1

    s =A@T
    s[0] = 0
    s[-1] = -2*dx
    u_1 = np.linalg.solve(B,s)
    return u_1
sol2 = Trapezoidal(T,D,x,dx,dt)

fig, ax = plt.subplots()
line, = ax.plot(x, sol2)
ax.set_title('Ecuación difusión 1D CrankNicholson')
ax.set_xlabel('x (m)')
ax.set_ylabel('T')



def animate(t):
    global T_desp
    T_desp = Trapezoidal(T_desp,D,x,dx,dt)
    line.set_ydata(T_desp)  
    return line,

ani = FuncAnimation(
    fig, animate, interval=60, blit=True)

plt.show()


#2D solution

D = 0.1
dx = x[1]-x[0]
x = np.linspace(-1, 1,50)
y = np.linspace(-1, 1,50)
X, Y = np.meshgrid(x, y)
T = 100 * np.exp(-20 * (X**2 + Y**2))
T_desp = np.copy(T)

alpha = D*dt/dx**2
def CrankNicholson2D(alpha,T):
    d = 1-2*alpha*np.ones(len(x))
    u = alpha*np.ones(len(x)-1)
    A = sparse.diags([d,u,u], [0,-1,1], shape = (len(x),len(x))).toarray()
    A[-1,-2] = -2*alpha
    A[-1,-1] = 1+2*alpha
    A[0,:] = A[-1,:] = 0 
    A[0,0] = A[-1,-1] = 1

    o = 1+2*alpha*np.ones(len(x))
    B = sparse.diags([o,-u,-u],[0,1,-1], shape = (len(x),len(x))).toarray()
    B[0,:] = B[-1,:] = 0
    B[0,0] = B[-1,-1] = 1
    
    I = sparse.eye(len(x)).toarray()
    bloquesA = [[None for _ in range(len(x))] for _ in range(len(x))]    
    bloquesB = [[None for _ in range(len(x))] for _ in range(len(x))]

    for i in range(0,len(x)):
        for j in range(0,len(x)):
            if i == j:    
                bloquesA[i][j] = A
                bloquesB[i][j] = B
            elif abs(i-j) == 1:
                bloquesA[i][j] = I*alpha
                bloquesB[i][j] = I*-alpha
            else:
                bloquesA[i][j] = sparse.csr_matrix((len(x),len(x)))
                bloquesB[i][j] = sparse.csr_matrix((len(x),len(x)))
            
    Abuena = sparse.bmat(bloquesA, format="csr").toarray()
    Bbuena =  sparse.bmat(bloquesB, format="csr").toarray()
    T_new = T.flatten()

    s = Abuena@T_new 
    r = np.linalg.solve(Bbuena,s)
    return r

sol = CrankNicholson2D(alpha,T)
rbuena = np.reshape(np.ravel(sol),(len(x),len(x)))

fig, ax = plt.subplots()
cax = ax.contourf(X,Y,rbuena, 20, cmap='inferno')
ax.set_title('Animación 2D cambiando valores de t = 0.001s')
ax.set_xlabel('x (m)')
ax.set_ylabel('T')
def update(frame):
    global T_desp, rbuena
    sol = CrankNicholson2D(alpha,T_desp)
    rbuena = np.reshape(np.ravel(sol),(len(x),len(x)))
    T_desp = np.reshape(sol, (len(x), len(x)))
    ax.clear()
    cax = ax.contourf(X, Y, rbuena, 20, cmap='inferno')
    return cax, 

anim = FuncAnimation(fig, update, frames=100, blit=False)
plt.show()
