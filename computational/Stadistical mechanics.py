# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 21:08:22 2024

@author: Usr

"""
#Stadistical mechanics for a crystal system bcc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


tipo = 'fcc'
try:
    coordx = int(input('Longitud de la caja de simulación en x (en parámetros de red) = '))
    coordy = int(input('Longitud de la caja de simulación en y (en parámetros de red) = '))
    coordz = int(input('Longitud de la caja de simulación en z (en parámetros de red) = '))
    condicion = int(input('Elige el tipo de condición: Periódica(1) o Libre(2), (elegir 1 o 2): '))
    

except:
    print('Error en la selección de uno o varios parámetros, inténtelo de nuevo')
    coordx = int(input('Longitud de la caja de simulación en x (en parámetros de red) = '))
    coordy = int(input('Longitud de la caja de simulación en y (en parámetros de red) = '))
    coordz = int(input('Longitud de la caja de simulación en z (en parámetros de red) = '))
    condicion = int(input('Elige el tipo de condición: Periódica(1) o Libre(2), (elegir 1 o 2): '))
    
dim = 3 #ya lo he puesto en 2 dimensiones para que sea LXL
pos = 0
a = 3.603
Lx = a*coordx
Ly = a*coordy
Lz = a*coordz
matrizcoord = np.zeros((coordx,coordy,coordz,dim))
#bcc
if tipo == 'fcc':
    cara1 = np.zeros((coordx,coordy,coordz,dim))
    cara2 = np.zeros((coordx,coordy,coordz,dim))
    cara3 = np.zeros((coordx,coordy,coordz,dim))
    for k in range(coordz):
        for j in range(coordy):
            for i in range(coordx):
              
                matrizcoord[i,j,k] = a*np.array([i,j,k])
               
                cara1[i,j,k] = a * np.array([i + 0.5, j + 0.5, k])
              
                cara2[i,j,k] = a * np.array([i, j + 0.5, k + 0.5])
              
                cara3[i,j,k] = a * np.array([i + 0.5, j, k + 0.5])
                
    matriz_total = np.concatenate((matrizcoord, cara1, cara2, cara3), axis=0)
    centro = np.mean(matriz_total, axis=(0, 1, 2))
   
    matrizcoord -= centro
    cara1 -= centro
    cara2 -= centro
    cara3 -= centro
    atomos_todos = np.concatenate((matrizcoord, cara1, cara2, cara3), axis=0)
   
    atomos_todos_bueno = atomos_todos.reshape(-1,3)
   
    atomos_x = atomos_todos_bueno[:,0]
    atomos_y = atomos_todos_bueno[:,1]
    atomos_z = atomos_todos_bueno[:,2]
    longitud = (len(atomos_todos_bueno))
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    ax.scatter(matrizcoord[:,:,:,0], matrizcoord[:,:,:,1], matrizcoord[:,:,:,2])
    ax.scatter(cara1[:,:,:,0], cara1[:,:,:,1], cara1[:,:,:,2]) 
    ax.scatter(cara2[:,:,:,0], cara2[:,:,:,1], cara2[:,:,:,2])
    ax.scatter(cara3[:,:,:,0], cara3[:,:,:,1], cara3[:,:,:,2])              
    ax.set_title('Red fcc')
    ax.set_xlabel('x (Ámstrongs)')
    ax.set_ylabel('y (Ámstrongs)')
    ax.set_zlabel('z (Ámstrongs)')
    plt.show()
"la lógica del programa es este: yo meto un cuadrado con los átomos, donde a cada lista le asigno un tipo \
de partícula, las de 1 cara, las de otra, etc, con el np.concatenate las junto todas en una matriz donde están \
todos los átomos (x,y) por las dimensiones de la caja, por tanto en el reshape se ponen los átomos en una matriz de 2 \
dimensionesxN átomos, no sé si se entiende muy bien pero esto no haría falta tocarlo porque ya está hecho para almacenar todos \
los átomos en una matriz de 2N"
           

def Energialibre(atomos_x,atomos_y,atomos_z, radiocorte,n,m,epsilon, longitud):
    
    potencial = []

    potencial_3d = np.zeros(longitud)
 
    fuerza = np.zeros(longitud)
    
    for i in range(len(atomos_x)):
        for j in range(i+1,len(atomos_x)):
            distancia = abs(np.sqrt((atomos_x[i]-atomos_x[j])**2+(atomos_y[i]-atomos_y[j])**2+(atomos_z[i]-atomos_z[j])**2))
       
            if distancia <= radiocorte:
                potencial_for = 4*epsilon*(((-sigma**n)/distancia**(n+1)) + (sigma**m)/(distancia**(m+1)))
                
                potencial.append((4*epsilon*((sigma/distancia)**n - (sigma/distancia)**m)))
                fuerza[i] += (-potencial_for*(atomos_x[i]-atomos_x[j])/distancia)
                fuerza[j] += (-potencial_for*(atomos_y[i]-atomos_y[j])/distancia)
              
                potencial_3d[i] += ((4*epsilon*((sigma/distancia)**n - (sigma/distancia)**m)))    
                potencial_3d[j] += ((4*epsilon*((sigma/distancia)**n - (sigma/distancia)**m)))   
    
    
    potencial = np.array(potencial)
 
    return potencial,potencial_3d,fuerza

#condiciones periodicas de contorno
def Periodicas(atomos_x,atomos_y,atomos_z,radiocorte,n,m,epsilon,Lx,Ly,Lz,longitud):

    potencial = []

    potencial_3d = np.zeros(longitud)
    fuerza = np.zeros(longitud)
    #aquí se hace lo de los vecinos, ya está tomado como se explicó en estadística con una suma primero hasta N y después i!= j hasta N
    #no debería dar problemas que vaya hasta la longitud de átomos_x porque ya es N átomos en cada dimensión, como va a ser cuadrada ya estaría bien
    
    for i in range(len(atomos_x)):
        for j in range(i+1,len(atomos_x)):
            deltax = atomos_x[j]-atomos_x[i]
            deltay = atomos_y[j]-atomos_y[i]
            deltaz = atomos_z[j]-atomos_z[i]
            if deltax > Lx/2:
                deltax = deltax - Lx
            if deltay > Ly/2:
                deltay = deltay - Ly
            if deltaz > Lz/2:
                deltaz = deltaz - Lz
            if deltax < -Lx/2:
                deltax = deltax + Lx
            if deltay < -Ly/2:
                deltay = deltay + Ly
            if deltaz < -Lz/2:
                deltaz = deltaz + Lz
            distancia = abs(np.sqrt(deltax**2+deltay**2+deltaz**2))
            
            if distancia <= radiocorte:
                potencial_for = 4*epsilon*(((-sigma**n)/distancia**(n+1)) + (sigma**m)/(distancia**(m+1)))
                potencial.append((4*epsilon*((sigma/distancia)**n - (sigma/distancia)**m)))
                fuerza[i] += (-potencial_for*(atomos_x[i]-atomos_x[j])/distancia)
                fuerza[j] += (-potencial_for*(atomos_y[i]-atomos_y[j])/distancia)
                
                potencial_3d[i] += ((4*epsilon*((sigma/distancia)**n - (sigma/distancia)**m)))    
                potencial_3d[j] += ((4*epsilon*((sigma/distancia)**n - (sigma/distancia)**m)))  

    potencial = np.array(potencial)
    
    return potencial,potencial_3d,fuerza

#datos
#aquí la temperatura debería entrar como dato pero no tengo claro qué hace falta
n = 12
m = 6
sigma = 2.3151
epsilon = 0.167
radiocorte = 3*sigma

if condicion == 1:
    Energia, potencial_3d,fuerza = Energialibre(atomos_x,atomos_y,atomos_z, radiocorte,n,m,epsilon,longitud)

    Energiatotlibre = np.sum(Energia)
    
    Energia_atomo = Energiatotlibre/len(atomos_todos_bueno)
    atomos_x = atomos_todos_bueno[:,0]
    atomos_y = atomos_todos_bueno[:,1]
    atomos_z = atomos_todos_bueno[:,2]
    color_u = np.array((potencial_3d-min(potencial_3d))/(max(potencial_3d)-min(potencial_3d)))
    color_f = np.array((fuerza-min(fuerza))/(max(fuerza)-min(fuerza)))
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    
    scatter = ax.scatter(atomos_x,atomos_y,atomos_z, c = color_u, cmap = 'viridis')
    cbar2 = plt.colorbar(scatter, ax=ax)
    cbar2.set_label('Magnitud del potencial')
    ax.set_title('Red fcc potencial (eV)')
    ax.set_xlabel('x (Ámstrongs)')
    ax.set_ylabel('y (Ámstrongs)')
    ax.set_zlabel('z (Ámstrongs)')
    plt.show()
    plt.close()
    
    fig = plt.figure()
    ax2 = fig.add_subplot(projection='3d')
    scatter = ax2.scatter(atomos_x,atomos_y,atomos_z, c = color_f, cmap = 'inferno')
    cbar2 = plt.colorbar(scatter, ax=ax2)
    cbar2.set_label('Magnitud de la fuerza')
    ax2.set_title('Red fcc fuerza (N)')
    ax2.set_xlabel('x (Ámstrongs)')
    ax2.set_ylabel('y (Ámstrongs)')
    ax2.set_zlabel('z (Ámstrongs)')
    plt.show()
    plt.close()
    
        
if condicion == 2:
        
    Energia, potencial_3d,fuerza = Periodicas(atomos_x,atomos_y,atomos_z,radiocorte,n,m,epsilon,Lx,Ly,Lz,longitud)
    
    Energiatotper = np.sum(Energia)
    
    
    Energia_atomo = Energiatotper/len(atomos_todos_bueno)
   
    atomos_x = atomos_todos_bueno[:,0]
    atomos_y = atomos_todos_bueno[:,1]
    atomos_z = atomos_todos_bueno[:,2]
    color_u = np.array((potencial_3d-min(potencial_3d))/(max(potencial_3d)-min(potencial_3d)))
    color_f = np.array((fuerza-min(fuerza))/(max(fuerza)-min(fuerza)))
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    
    scatter = ax.scatter(atomos_x,atomos_y,atomos_z, c = color_u, cmap = 'viridis')
    cbar2 = plt.colorbar(scatter, ax=ax)
    cbar2.set_label('Magnitud del potencial')
    ax.set_title('Red fcc potencial (eV)')
    ax.set_xlabel('x (Ámstrongs)')
    ax.set_ylabel('y (Ámstrongs)')
    ax.set_zlabel('z (Ámstrongs)')
    plt.show()
    plt.close()
    
    fig = plt.figure()
    ax2 = fig.add_subplot(projection='3d')
    scatter = ax2.scatter(atomos_x,atomos_y,atomos_z, c = color_f, cmap = 'inferno')
    cbar2 = plt.colorbar(scatter, ax=ax2)
    cbar2.set_label('Magnitud de la fuerza')
    ax2.set_title('Red fcc fuerza (eV/A)')
    ax2.set_xlabel('x (Ámstrongs)')
    ax2.set_ylabel('y (Ámstrongs)')
    ax2.set_zlabel('z (Ámstrongs)')
    plt.show()
    plt.close()

#EXTRA 1
#aquí se trabaja que la temperatura varía en función de la energía cinética media del sistema
T = float(input('Introduce la temperatura (en K): '))

matrizcoord = np.random.random((coordx,coordy,coordz,dim)) -0.5
cara1 = np.random.random((coordx,coordy,coordz,dim)) -0.5
cara2 = np.random.random((coordx,coordy,coordz,dim)) -0.5
cara3 = np.random.random((coordx,coordy,coordz,dim)) -0.5
atomos_todos = np.concatenate((matrizcoord, cara1, cara2, cara3), axis=0)
atomos_todos_bueno = atomos_todos.reshape(-1,3)
atomos_x = atomos_todos_bueno[:,0]
atomos_y = atomos_todos_bueno[:,1]
atomos_z = atomos_todos_bueno[:,2]

vx = (np.random.random(len(atomos_x)) -0.5)*1e-10
vy = (np.random.random(len(atomos_y)) - 0.5)*1e-10
vz = (np.random.random(len(atomos_z)) -0.5)*1e-10
N = len(atomos_x)+len(atomos_y)+len(atomos_z)
m_cobre = 63.55
Ec = (1/2)*m_cobre*(vx**2+vy**2+vz**2)
Kb = 8.611024e-5*1.602177e-19
Trandom = (2/3)*np.mean(Ec)/(Kb)
sc = np.sqrt(T/Trandom)
vx_new = vx*sc
vy_new = vy*sc
vz_new = vz*sc

Ec = (1/2)*m_cobre*(vx_new**2+vy_new**2+vz_new**2)
Trandom2 = (2/3)*np.mean(Ec)/(Kb)
#obtenemos un vector con la temperatura que esperábamos, imprimimos el valor medio
print('La temperatura pedida para este programa es: ',np.mean(Trandom2))

#EXTRA 2
#esto seguramente no hace falta pero lo dejo por si acaso se puede aprovechar algo
#no lo tengo bien lo de sacar las fuerzas 
nt = int(input('Nº de pasos de la simulación: '))
k = 1
h = 1e-15
m_kg = m_cobre*1.66054e-27
Kb_eV = 8.611024e-5
Ec = []
Ep = []
Temp = []

def f(r):
    
    a = -k*r
    return a

r = np.random.uniform(-5,5,(len(atomos_x),3))*1e-10


v = np.zeros((len(vx_new),3))

v[:,0] = vx_new
v[:,1] = vy_new
v[:,2] = vz_new

a = f(r)/m_kg
Ep0 = (0.5*k*np.linalg.norm(r, axis = 1)**2)/1.602177e-19
Epi = np.sum(Ep0)
Ec0 = (0.5*m_kg*np.linalg.norm(v,axis = 1)**2)/1.602177e-19 
Temp0 = (2/3)*np.mean(Ec0)/(Kb_eV)
Ep.append(np.sum(Ep0)-Epi)
Ec.append(np.sum(Ec0))
Temp.append(Temp0)

for i in range(0,nt):  
    r_2 = r + (h*v) + 0.5*a*h**2
    a_2 = f(r_2)/m_kg
    v_2 = v + 0.5*h*(a+a_2)

    Epfor = (0.5*k*np.linalg.norm(r_2, axis = 1)**2)/1.602177e-19
    Ecfor = (0.5*m_kg*np.linalg.norm(v_2,axis = 1)**2)/1.602177e-19  
    Tempfor = (2/3)*np.mean(Ecfor)/(Kb_eV)
    Ep.append(np.sum(Epfor)-Epi)
    Ec.append(np.sum(Ecfor))
    Temp.append(Tempfor)
    r,v,a = r_2,v_2,a_2


x = np.linspace(0,nt+1,nt+1) 
Ec = np.array(Ec)
Ep = np.array(Ep)
Temp = np.array(Temp)
plt.title('Energía total para una partícula')
plt.plot(x,Ec, label = 'Energía cinética')
plt.plot(x,Ep, label = 'Energía potencial')
plt.plot(x,Ec+Ep, label = 'Energía total')
plt.xlabel('número de pasos')
plt.ylabel('E (eV)')
plt.legend(loc = 'best')
plt.show()

plt.plot(x,Temp)
plt.title('Evolución de la temperatura para una partícula')
plt.xlabel('Nº de pasos')
plt.ylabel('T(K)')

plt.show()

print('La presentación de la gráfica, tanto la de la energía como la temperatura nos indica que los átomos están vibrando dentro de la celda.')
