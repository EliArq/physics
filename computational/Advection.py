
#Centered differences for advection equation and burgers equation
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

x = np.linspace(0,5,100)
c = 2
u = np.zeros(100)
dx = x[1] -x[0]
nt = 1000
dt = 0.01

for j in range(0,len(x)-1):
    if 1 <= x[j] <= 2:
        u[j] = 1
        
def DiferenciasCentradas(c,u,dx,dt):
    factor = (c*dt)/(2*dx)
    u_2 = np.copy(u)
    for j in range(1,len(u)-1):
        u_2[j] = u[j] - factor*(u[j+1]-u[j-1])
    u_2[0] = u[0]-factor*(u[1]-u[-2])
    u_2[-1] = u_2[0]     
    return u_2

def Downwind(c,u,dx,dt):
    factor = (c*dt)/dx
    u_3 = np.copy(u)
    for j in range(1,len(u)-1):
        u_3[j] = u[j] - factor*(u[j+1]-u[j])
    u_3[0] = u[0] - factor*(u[1]-u[0])
    u_3[-1] = u_3[0]
    
    return u_3

def Upwind(c,u,dx,dt):
    factor = (c*dt)/dx
    u_4 = np.copy(u)
    for j in range(1,len(u)):
        u_4[j] = u[j] - factor*(u[j]-u[j-1])
    
    u_4[0] = u_4[-1]
    return u_4

u_centradas, u_downwind, u_upwind = np.copy(u), np.copy(u), np.copy(u)

fig, axs = plt.subplots(3,1,figsize = (8,10))
titles = ['Diferencias centradas','Downwind','Upwind']
lines = []
for ax,title in zip(axs,titles):
    line, = ax.plot(x, u)
    ax.set_title(title)
    ax.set_xlabel('x (m)')
    ax.set_ylabel('u(x)')
    ax.set_ylim(-0.5,1.5)
    lines.append(line)

def animate(t):
    global u_centradas, u_downwind, u_upwind
    u_centradas = DiferenciasCentradas(c,u_centradas,dx,dt)
    u_downwind = Downwind(c, u_downwind, dx,dt)
    u_upwind = Upwind(c,u_upwind,dx,dt)
    lines[0].set_ydata(u_centradas)
    lines[1].set_ydata(u_downwind)
    lines[2].set_ydata(u_upwind)

    return lines

ani = FuncAnimation(
    fig, animate, frames = 1000, interval=30, blit=True)

plt.tight_layout()
plt.show()
plt.close()

x1 = np.linspace(0,1,100)
dx1 = x1[1]-x1[0]

dt1 = 0.002
nt = 10000
u1= 2 + 0.5*np.sin(2*np.pi*x1)
#if dt < (dx/np.max(u)):
 #   print('cumple condicion')
    
def BurgersUpwind(dx1,dt1,u1):
    u_5 = np.copy(u1)
    for i in range(1,len(u1)):
        u_5[i] = u1[i] - (0.5*(u1[i]**2)-0.5*(u1[i-1]**2))*dt1/dx1
    u_5[0] = u_5[-1]      
    return u_5


def BurgersCentradas(dx1,dt1,u1):
    u_3 = np.copy(u1)
    for j in range(1,len(u1)-1):
        u_3[j] = u_3[j] - u1[j]*(u1[j+1]-u1[j-1])*dt1/(2*dx1)
    u_3[-1] = u_3[0]
        
    return u_3

u_centradas_bur, u_upwind_bur = np.copy(u1), np.copy(u1)

fig, axs = plt.subplots(2,1,figsize = (8,10))
titles = ['Burgers Diferencias centradas','Burgers Upwind']
lines = []
for ax,title in zip(axs,titles):
    line, = ax.plot(x1, u1)
    ax.set_title(title)
    ax.set_xlabel('x (m)')
    ax.set_ylabel('u(x,t)')
    lines.append(line)

def animate(t):
    global u_centradas_bur, u_upwind_bur
    u_centradas_bur = BurgersCentradas(dx1,dt1,u_centradas_bur)
    u_upwind_bur = BurgersUpwind(dx1,dt1,u_upwind_bur)
    lines[0].set_ydata(u_centradas_bur)
    lines[1].set_ydata(u_upwind_bur)

    return lines
ani = FuncAnimation(
    fig, animate, frames = 1000, interval=30, blit=True)

plt.tight_layout()
plt.show()
plt.close()

t = np.linspace(0,1,100)

def Energia(u1):
    E = 0.5*sum(u1**2)*dx1
    return E

E = np.zeros(len(t))
for i in range(len(t)):
    u1 = BurgersUpwind(dx1,dt1,u1)
    E[i] = Energia(u1)


plt.plot(t,E)
plt.title('Energía en función del tiempo Burgers Upwind (estable)')
plt.xlabel('t(s)')
plt.ylabel('E(J)')
plt.show()
plt.close()
print('Como se puede observar, el método es estable pero presenta una gran cantidad de dispersión, debida al gran incremento de tiempo dt = 0.002')


