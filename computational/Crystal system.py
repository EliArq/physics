# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#Crystal system construction: cubic, fcc, bcc and diamond
#L = 1
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

tipo = input('Tipo de red(bcc,fcc,cubic,diamond): ')
coordx = int(input('Longitud de la caja de simulación en x (en parámetros de red) = '))
coordy = int(input('Longitud de la caja de simulación en y (en parámetros de red) = '))
coordz = int(input('Longitud de la caja de simulación en z (en parámetros de red) = '))
a = int(input('Dimensiones de la red = '))
dim = 3
matrizcoord = np.zeros((coordx,coordy,coordz,dim))
pos = 0
file_name = 'Archivo salida.txt'
file = open(file_name, 'w')
#cubica
if tipo == 'cubic':
    
    for k in range(0,coordz):
        for j in range(0,coordy):
            for i in range(0,coordx):
                matrizcoord[i,j,k,:] = a*np.array([i,j,k])
                pos += 1
                file.write('f {pos},{matrizcoord[i,j,k]} \n')
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    ax.scatter(matrizcoord[:,:,:,0], matrizcoord[:,:,:,1], matrizcoord[:,:,:,2]) 
    ax.set_title('Red cúbica')
    plt.show()


#bcc
if tipo == 'bcc':
   
    medio = np.zeros((coordx,coordy,coordz,dim))
    for k in range(coordz):
        for j in range(coordy):
            for i in range(coordx):
                matrizcoord[i,j,k] = a*np.array([i,j,k])
                medio[i,j,k] = a * np.array([i + 0.5, j + 0.5, k + 0.5])
                
                
    matriz_total = np.concatenate((matrizcoord, medio), axis=0)
    centro = np.mean(matriz_total, axis=(0, 1, 2))
    matrizcoord -= centro
    medio -= centro
    
    for k in range(coordz):
        for j in range(coordy):
            for i in range(coordx):
                pos += 2
                file.write(f'{pos},{matrizcoord[i,j,k]},{medio[i,j,k]} \n')
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    ax.scatter(matrizcoord[:,:,:,0], matrizcoord[:,:,:,1], matrizcoord[:,:,:,2])
    ax.scatter(medio[:,:,:,0], medio[:,:,:,1], medio[:,:,:,2]) 
    ax.set_title('Red bcc')
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')           
    plt.show()
                

#fcc
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
    
    for k in range(coordz):
        for j in range(coordy):
            for i in range(coordx):
                pos += 4
                file.write(f'{pos}, {matrizcoord[i,j,k]},{cara1[i,j,k]},{cara2[i,j,k]},{cara3[i,j,k]} \n')
                
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
                
                
#diamante
if tipo == 'diamond':
    
    cara1, cara2, cara3, cara4, cara5, cara6, cara7 = np.zeros((coordx,coordy,coordz,dim)), np.zeros((coordx,coordy,coordz,dim)), np.zeros((coordx,coordy,coordz,dim)), np.zeros((coordx,coordy,coordz,dim)), np.zeros((coordx,coordy,coordz,dim)), np.zeros((coordx,coordy,coordz,dim)), np.zeros((coordx,coordy,coordz,dim))
    
    for k in range(coordz):
        for j in range(coordy):
            for i in range(coordx):
                matrizcoord[i,j,k] = a*np.array([i,j,k])
                cara1[i,j,k] = a * np.array([i + 0.5, j + 0.5, k])
                cara2[i,j,k] = a * np.array([i, j + 0.5, k + 0.5])
                cara3[i,j,k] = a * np.array([i + 0.5, j, k + 0.5])
                cara4[i,j,k] = a * np.array([i+0.25, j+0.25, k+0.25])
                cara5[i,j,k] = a * np.array([i+0.75, j+0.75, k+0.25])
                cara6[i,j,k] = a * np.array([i+0.75, j+0.25, k+0.75])
                cara7[i,j,k] = a * np.array([i+0.25, j+0.75, k+0.75])
                
    
    
    matriz_total = np.concatenate((matrizcoord, cara1, cara2, cara3, cara4, cara5, cara6, cara7), axis=0)
    centro = np.mean(matriz_total, axis=(0, 1, 2))


    matrizcoord -= centro
    cara1 -= centro
    cara2 -= centro
    cara3 -= centro
    cara4 -= centro
    cara5 -= centro
    cara6 -= centro
    cara7 -= centro
        
    for k in range(coordz):
        for j in range(coordy):
            for i in range(coordx):
                pos += 8
                file.write(f'{pos}, {matrizcoord[i,j,k]},{cara1[i,j,k]},{cara2[i,j,k]},{cara3[i,j,k]},{cara4[i,j,k]},{cara5[i,j,k]},{cara6[i,j,k]},{cara7[i,j,k]} \n')
    
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    ax.scatter(matrizcoord[:,:,:,0], matrizcoord[:,:,:,1], matrizcoord[:,:,:,2])
    ax.scatter(cara1[:,:,:,0], cara1[:,:,:,1], cara1[:,:,:,2]) 
    ax.scatter(cara2[:,:,:,0], cara2[:,:,:,1], cara2[:,:,:,2])
    ax.scatter(cara3[:,:,:,0], cara3[:,:,:,1], cara3[:,:,:,2]) 
    ax.scatter(cara4[:,:,:,0], cara4[:,:,:,1], cara4[:,:,:,2])  
    ax.scatter(cara5[:,:,:,0], cara5[:,:,:,1], cara5[:,:,:,2])  
    ax.scatter(cara6[:,:,:,0], cara6[:,:,:,1], cara6[:,:,:,2])  
    ax.scatter(cara7[:,:,:,0], cara7[:,:,:,1], cara7[:,:,:,2])               
    ax.set_title('Red diamond')
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    plt.show()

file.write(f'El número total de posiciones es: {pos} \n')
file.write(f'Las dimensiones de la red son: {coordx} , {coordy}, {coordz} \n')
file.close()
print('Archivo cerrado')

