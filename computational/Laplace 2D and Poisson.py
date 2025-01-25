# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 19:27:19 2024

@author: Usr
"""


import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
from scipy.sparse.linalg import spsolve

def TridiagonalSolver(d,o,u,r):
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
 
    return soluciones

#Poisson equation solve
N = 100
dominio = np.linspace(0,np.pi,N)
def Poisson(N):
    dominio = np.linspace(0,np.pi,N)    
    cond_iniciales = np.array([0,-np.pi**2])
    delta_x = (dominio[-1]-dominio[0])/(N-1)
    u = 1*np.ones(N-1)
    o = 1*np.ones(N-1)
    d = -2*np.ones(N)
    o[0] = 0
    d[0] = 1
    u[-1] = 0
    d[-1] = 1
    r = ((2-(dominio**2))*np.sin(dominio) + 4*dominio*np.cos(dominio))*delta_x**2
    r[0] = cond_iniciales[0]
    r[-1] -= (-delta_x * cond_iniciales[1])
    resultado = TridiagonalSolver(d, o, u, r)
    
    return resultado

y = Poisson(N)
plt.plot(dominio, y, label="Solución f(x)", color='blue')
plt.title("Solución de la ecuación de Poisson")
plt.xlabel("x")
plt.ylabel("f(x)")
plt.show()
plt.close()
        
#Laplace 2D solve using block matrix
L = 1
N = 50

d = -4*np.ones(N)
u = 1*np.ones(N)
o = 1*np.ones(N)
diag= sparse.diags([d,u,o], [0,-1,1], shape = (N,N)).toarray() 

identidad = sparse.eye(N).toarray()

def Laplace(L,N):
    x = np.linspace(0,L,N)
    r = np.zeros(N*N)
    for i in range(1,N+1):
        r[-N+i-1] = x[i-1]*(1-x[i-1])
   
    return r

r = Laplace(L,N)


derivadas = sparse.kron(identidad,diag).toarray() 
bloques = [[None for _ in range(N)] for _ in range(N)]


for i in range(N):
    bloques[i][i] = diag
    

for i in range(N - 1):
    bloques[i][i + 1] = sparse.eye(N)
    bloques[i + 1][i] = sparse.eye(N)

#bmat method for join block matrix
A = sparse.bmat(bloques, format="csr").toarray()


phi = spsolve(A,r)

phibuena = np.reshape(np.ravel(phi),(N,N))


plt.contourf(phibuena, 8, cmap=plt.cm.hot) 
C = plt.contour(phibuena, 8, colors='black') 
plt.clabel(C, inline=True) 

plt.xticks(())
plt.yticks(())


plt.figure()
plt.imshow(phibuena,cmap = 'inferno', origin = 'lower')
plt.colorbar()
plt.show()


#Jacobi method iterative
error = 1e-6
L = 1.
N = 50   
            
def Jacobi(error,N,L):
    x = np.linspace(0, L, N)

    phi = np.zeros((N, N))
    phi2 = np.zeros((N, N)) 
    phi[0, :] = 0             
    phi[-1, :] = 0            
    phi[:, 0] = 0 
    phi[:, -1] = x*(1-x)           
    phi2 = np.copy(phi)
    convergencia = 1
    while convergencia > error :
        phi2[1:-1,1:-1] = (1/4)*(phi[2:,1:-1] + phi[:-2,1:-1] + phi[1:-1,2:] + phi[1:-1,:-2])
        convergencia = np.max(np.absolute(phi2-phi))
        phi,phi2 = phi2,phi

    return phi


phi = Jacobi(error,N,L)

plt.contourf(phi, 10, alpha=.75, cmap=plt.cm.inferno)
C = plt.contour(phi, 10, colors='black') 
plt.clabel(C, inline=True, fontsize=20) 
plt.xticks(())
plt.yticks(())
plt.figure()
plt.imshow(phi, origin = 'lower')
plt.show()

    
        