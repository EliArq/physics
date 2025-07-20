# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#Metropolis MonteCarlo mechanical stadistics
import numpy as np
import matplotlib.pyplot as plt

N = int(input('Número de partículas: '))
beta = float(input('Valor de KT: '))
pasos = int(input('Número de pasos Monte Carlo: '))


def Metropolis(pasos,N,beta):
    L = 1
    m = 1
    pi = 1
    hbarra = 1
    nx = np.ones(N)
    ny = np.ones(N)
    nz = np.ones(N)
    E_1 = ((pi**2*hbarra**2*(nx**2+ny**2+nz**2))/(2*m*L**2))
    E_old = E_1
    energias = []
    #hay que comparar la E1 con un número que va entre 0 y 1
    for i in range(pasos):  
        pe = np.random.choice(N)
        direccion = np.random.choice(['x','y','z'])
        signo = np.random.choice([-1,1])
        R = np.random.rand()
        if direccion == 'x' and nx[pe] >= 1:
            nx[pe] += signo
        elif direccion == 'y' and ny[pe] >= 1:
            ny[pe] += signo      
        elif direccion == 'z' and nz[pe] >= 1:
            nz[pe] += signo
        E_1[pe] = ((pi**2 * hbarra**2 * (nx[pe]**2 + ny[pe]**2 + nz[pe]**2)) / (2 * m * L**2))    
        deltaE = E_1[pe]-E_old[pe]
        if deltaE < 0 or R < np.exp(-deltaE/(beta)):
            E_old[pe] += deltaE
        else:
            if direccion == 'x':
              nx[pe] -= signo
            elif direccion == 'y':
              ny[pe] -= signo
            elif direccion == 'z':
              nz[pe] -= signo
        energias.append(np.sum(E_old))

    return energias


listaE = Metropolis(pasos,N,beta)
listaE = np.array(listaE)
x = np.linspace(0,pasos,len(listaE))
totalE = np.sum(listaE)/len(listaE)
print(f'La energía promedio para un número {N} de partículas es {totalE:.3f} [E]')

plt.title('Energía vs Pasos Monte Carlo')
plt.plot(x,listaE)
plt.show()

#comparamos con varios valores de T y N
beta2 = np.linspace(beta,100,10)
compara = []
for j in range(len(beta2)):
    E = Metropolis(pasos,N,beta2[j])
    compara.append(np.mean(E))
    
plt.title('Energía total para diferente T')
plt.plot(beta2,compara, 'o-')
plt.show()


N2 = np.arange(N,10000,1000)
comparaN = []
for m in range(len(N2)):
    nu = Metropolis(pasos,N2[m],beta)
    comparaN.append(np.mean(nu))
        
plt.title('Energía en función del tamaño N')
plt.plot(N2,comparaN, 'o-')
plt.show()

print('En este caso, aunque no se aprecie porque no es preciso en la gráfica, \n'
      'la energía disminuye ligeramente al aumentar beta, ya que se restringe la accesibilidad de estados \n'
      'y, respecto al número de pasos, la energía aumenta al aumentar la cantidad de posibilidades de interacción \n'
      'existen dentro del sistema y la energía promedio aumenta')