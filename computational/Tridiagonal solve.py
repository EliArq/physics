# -*- coding: utf-8 -*
#Tridiagonal matrix using Thomas Algorithm
import numpy as np
import time
import matplotlib.pyplot as plt


def TridiagonalSolver(d,o,u,r):
    inicio = time.time()
    ro = np.zeros(len(r))
    n = len(d)
    ro[0] = r[0]/d[0]
    
    h = np.zeros(len(d)-1)
    h[0] = (o[0]/d[0])
    i = 1
    j = 1
    while i < n-1:
        h[i] = o[i]/(d[i]-(u[i]*h[i-1]))
        i+=1
    while j < n:
        ro[j] = (r[j]-(u[j-1]*ro[j-1]))/(d[j]-(u[j-1]*h[j-1]))
        j += 1

    soluciones = np.zeros(len(d))
    soluciones[-1] = ro[-1] 
    
    for i in range(n-2,-1,-1):
        soluciones[i] = ro[i] - (h[i]*soluciones[i+1])
    termino = time.time()
    return soluciones, termino,inicio
    

#la diagonal debe ser mayor a la suma de u y o para cada i
matriz = [[3,1,0,0],[1,5,2,0],[0,2,7,3],[0,0,3,9]]
d = np.diag(matriz)
u = np.diag(matriz,-1)

o = np.diag(matriz,1)
r = np.array([1,2,3,4])

inicio,termino = TridiagonalSolver(d,o,u,r)[1:]    
tiempo = termino-inicio 

#print(TridiagonalSolver(d,o,u,r))


#Time efficiency
longitud_matrix = []
tiempo_ejec = []

for i in range(1,10,3):
    n = 1000*i
    d = np.random.rand(n)
    u = np.random.rand(n-1)
    o = np.random.rand(n-1)
    matriz2 = np.diag(d) + np.diag(u,1) + np.diag(o,-1)
    r = np.random.rand(1000*(i))

    solucion2,termino,inicio = TridiagonalSolver(d,o,u,r)
    tiempo = termino-inicio
    tiempo_ejec.append(tiempo)
    longitud_matrix.append(len(matriz2))


plt.plot(longitud_matrix,tiempo_ejec, 'ro')
plt.ylabel('tiempo de ejecución (s)')
plt.xlabel('longitud de la matriz (len)')
plt.title('Representación del tiempo de ejecución respecto del tamaño de la matriz')
plt.show()


#Linalg comparison 
longitud_matrix2 =[]
tiempo_ejec2 = []
longitud_matrix3 = []
tiempo_ejec3 = []

for m in range(1,10,3):
    n = 1000*m
    d = np.random.rand(n)
    u = np.random.rand(n-1)
    o = np.random.rand(n-1)
    matriz3 = np.diag(d) + np.diag(u,1) + np.diag(o,-1)
    r = np.random.rand(1000*(m))
    inicio = time.time()

    matriz_x = np.linalg.solve(matriz3,r)
    termino = time.time()    
    tiempo = termino-inicio 
    tiempo_ejec2.append(tiempo)
    longitud_matrix2.append(len(matriz_x))
    
for j in range(1,10,3):
    n = 1000*j
    d = np.random.rand(n)
    u = np.random.rand(n-1)
    o = np.random.rand(n-1)
    matriz3 = np.diag(d) + np.diag(u,1) + np.diag(o,-1)
    r = np.random.rand(1000*(j))
    inicio = time.time()

    matriz3_inv = np.linalg.inv(matriz3)

    sol1 = matriz3_inv@r
    termino = time.time()    
    tiempo = termino-inicio 

    tiempo_ejec3.append(tiempo)
    longitud_matrix3.append(len(matriz3))


plt.plot(tiempo_ejec, longitud_matrix, 'ro', 
         tiempo_ejec2, longitud_matrix2, 'yo', 
         tiempo_ejec3, longitud_matrix3, 'bo')


plt.legend(labels=['tridiagonal', 'linalg', 'inversa'], loc = 'best')
plt.xlabel('tiempo de ejecución (s)')
plt.ylabel('longitud de la matriz (len)')
plt.title('Representación del tiempo de ejecución respecto del tamaño de la matriz')
plt.show()
#print('Teniendo en cuenta que el color rojo corresponde al método "a mano", mientras que el amarillo corresponde a numpy solve y el azul a multiplicar por la inversa, el método más ineficiente es multiplicar la inversa de la matriz, y el más eficiente es el hecho "a mano"')