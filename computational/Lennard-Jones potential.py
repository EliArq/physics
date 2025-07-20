# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 16:48:15 2024

@author: alumno
"""

#Lennard-Jones potential 
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
    
dim = 3
pos = 0
a = 3.603
Lx = a*coordx
Ly = a*coordy
Lz = a*coordz
matrizcoord = np.zeros((coordx,coordy,coordz,dim))
#bcc
if tipo == 'fcc':
    atomos_centro = []
    atomos_lado1 = []
    atomos_lado2 = []
    atomos_lado3 = []
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
    #hemos juntado todos los átomos en celdas donde tenemos 8 átomos por celda, 2x2 celdas unitarias y las 3 dimensiones
    #tenemos que redimensionar la matriz para que todos los átomos tengan sus coordenadas y su dimensión correspondiente
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
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    plt.show()
                
#Boundary conditions
def Energialibre(atomos_x,atomos_y,atomos_z, radiocorte,n,m,epsilon, longitud):
    
    potencial = []

    potencial_3d = np.zeros(longitud)
    
    for i in range(len(atomos_x)):
        for j in range(i+1,len(atomos_x)):
            distancia = abs(np.sqrt((atomos_x[i]-atomos_x[j])**2+(atomos_y[i]-atomos_y[j])**2+(atomos_z[i]-atomos_z[j])**2))

            if distancia <= radiocorte:
                potencial.append((4*epsilon*((sigma/distancia)**n - (sigma/distancia)**m)))
                potencial_3d[i] += ((4*epsilon*((sigma/distancia)**n - (sigma/distancia)**m)))    
                potencial_3d[j] += ((4*epsilon*((sigma/distancia)**n - (sigma/distancia)**m)))   
    
    potencial = np.array(potencial)

    return potencial,potencial_3d

def Periodicas(atomos_x,atomos_y,atomos_z,radiocorte,n,m,epsilon,Lx,Ly,Lz,longitud):

    potencial = []

    potencial_3d = np.zeros(longitud)
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
                potencial.append((4*epsilon*((sigma/distancia)**n - (sigma/distancia)**m)))
                potencial_3d[i] += ((4*epsilon*((sigma/distancia)**n - (sigma/distancia)**m)))    
                potencial_3d[j] += ((4*epsilon*((sigma/distancia)**n - (sigma/distancia)**m)))  
                
    
    potencial = np.array(potencial)
    
    return potencial,potencial_3d

#datos
n = 12
m = 6
sigma = 2.3151
epsilon = 0.167
radiocorte = 3*sigma

if condicion == 1:
    Energia = Energialibre(atomos_x,atomos_y,atomos_z, radiocorte,n,m,epsilon,longitud)[0]
    potencial_3d = Energialibre(atomos_x,atomos_y,atomos_z, radiocorte,n,m,epsilon,longitud)[1]

    Energiatotlibre = np.sum(Energia)
    print(f'La energía total del sistema es (Ev): {Energiatotlibre:.3f}')
    Energia_atomo = Energiatotlibre/len(atomos_todos_bueno)
    print(f'La energía por unidad de átomo es (Ev/átomo): {Energia_atomo:.3f}')
    print('Para esta condición de contorno, tenemos una energía menor ya que consideramos que la red tiene una longitud \n'
          f'{longitud}x{dim}, por tanto la aportación de los átomos termina en los límites de la red.')
    atomos_x = atomos_todos_bueno[:,0]
    atomos_y = atomos_todos_bueno[:,1]
    atomos_z = atomos_todos_bueno[:,2]
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(atomos_x,atomos_y,atomos_z, c = potencial_3d, cmap = 'inferno')
    ax.set_title('Red fcc potencial (eV)')
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    plt.show()
    
        
if condicion == 2:
        
    Energia = Periodicas(atomos_x,atomos_y,atomos_z,radiocorte,n,m,epsilon,Lx,Ly,Lz,longitud)[0]
    potencial_3d = Periodicas(atomos_x,atomos_y,atomos_z,radiocorte,n,m,epsilon,Lx,Ly,Lz,longitud)[1]

    Energiatotper = np.sum(Energia)
    print(f'La energía total del sistema es (Ev): {Energiatotper:.3f}')
    Energia_atomo = Energiatotper/len(atomos_todos_bueno)
    print(f'La energía por unidad de átomo es (Ev/átomo): {Energia_atomo:.3f}')
    print('Para esta condición de contorno, tenemos una energía mayor ya que consideramos que la red es "infinita" \n'
          'por tanto existe mayor interacción entre los átomos, tanto dentro como fuera de las celdas, generando así mayor potencial.')
    atomos_x = atomos_todos_bueno[:,0]
    atomos_y = atomos_todos_bueno[:,1]
    atomos_z = atomos_todos_bueno[:,2]
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(atomos_x,atomos_y,atomos_z, c = potencial_3d, cmap = 'inferno')
    ax.set_title('Red fcc potencial (eV)')
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    plt.show()
    
#modificamos el radio de corte y comparamos los resultados

#tomamos un array de valores del radio de corte y comparamos los resultados
radiocorte = np.array([sigma/2,sigma,5*sigma,10*sigma])
for i in radiocorte:
        Energia = Energialibre(atomos_x,atomos_y,atomos_z, i,n,m,epsilon,longitud)[0]
        potencial_3d = Energialibre(atomos_x,atomos_y,atomos_z, i,n,m,epsilon,longitud)[1]
        Energiatotlibre = np.sum(Energia)
        print(f'La energía total del sistema para el radio de corte {i:.3f} (A) es (Ev) en condición libre: {Energiatotlibre:.3f}')
        Energia_atomo = Energiatotlibre/len(atomos_todos_bueno)
        print(f'La energía por unidad de átomo es (Ev/átomo): {Energia_atomo:.3f}')
       
        Energia = Periodicas(atomos_x,atomos_y,atomos_z,i,n,m,epsilon,Lx,Ly,Lz,longitud)[0]
        potencial_3d = Periodicas(atomos_x,atomos_y,atomos_z,i,n,m,epsilon,Lx,Ly,Lz,longitud)[1]
        Energiatotper = np.sum(Energia)
        print(f'La energía total del sistema para el radio de corte {i:.3f} (A) es (Ev) en condición periódica: {Energiatotper:.3f}')
        Energia_atomo = Energiatotper/len(atomos_todos_bueno)
        print(f'La energía por unidad de átomo es (Ev/átomo): {Energia_atomo:.3f}')

print('En este caso podemos observar que existe un límite donde el potencial es 0 pues la distancia del radio de corte es tan pequeña que no es  \n'
      'posible obtener ningún átomo cerca de otro que puedan generar ningún tipo de potencial, independientemente del tipo de condición que se imponga a la red.')
print('Para los casos donde se toman distancias mayores tenemos resultados similares, ya que al final tomamos distancias tan grandes que estamos tomando \n'
      'muchos átomos, y, bajo condiciones periódicas, la energía total puede aproximarse a una contribución promedio por átomo multiplicada por el número total de átomos \n'
      'considerados en la celda.')