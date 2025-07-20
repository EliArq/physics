# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 17:10:44 2024

@author: alumno
"""
#MonteCarlo integrations method 
import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return (np.sin(1/(x*(2-x))))**2

valormax = float(input('Introduce el valor máximo de la función: '))
nt = int(input('Introduce el número de puntos de Monte Carlo: '))

x = np.linspace(0,2,nt)
a = x[0]
b = x[-1]
plt.plot(x, f(x))
plt.title('Valor de la función entre $x = 0$ y $x = 2$')
plt.xlabel('Número de puntos')
plt.ylabel('f(x)')
plt.show()
plt.close()

def MonteCarlo(valormax,nt,a,b):
    y = np.random.uniform(0,valormax,nt)
    x = np.random.uniform(a,b,nt) 

    x_dentro = []
    y_dentro = []
    x_fuera = []
    y_fuera = []
    for i in range(0,nt):
        valor = f(x[i])
        if y[i] < valor:
            x_dentro.append(x[i])
            y_dentro.append(y[i])
        else:
            x_fuera.append(x[i])
            y_fuera.append(y[i])
    
    x_dentro,y_dentro = np.array(x_dentro),np.array(y_dentro)
    frac_pares = len(x_dentro)/nt
    I = valormax*(b-a)*frac_pares       
    desv = np.sqrt(I*(valormax*(b-a)-I))/np.sqrt(nt)
    return I,desv,x_dentro,y_dentro,x_fuera,y_fuera

resultados = MonteCarlo(valormax,nt,a,b)
I,desv,x_dentro,y_dentro,x_fuera,y_fuera = resultados
print(f'El valor de la integral es {I:.3f} ± {desv:.3f}')


#3
valores_I = []
n_pasos = []
desv_valores = []
#a partir de 100k tarda demasiado en hacer los cálculos
if nt <= 5000:
    for i in range(1,nt,10):
        n_pasos.append(i)
        I,desv = MonteCarlo(valormax,i,a,b)[0], MonteCarlo(valormax,i,a,b)[1]
        valores_I.append(I)
        desv_valores.append(desv)
if nt >= 100000:
    for i in range(1,nt,10000):
        n_pasos.append(i)
        I,desv = MonteCarlo(valormax,i,a,b)[0], MonteCarlo(valormax,i,a,b)[1]
        valores_I.append(I)
        desv_valores.append(desv)
else:
    for i in range(1,nt,1000):
        n_pasos.append(i)
        I,desv = MonteCarlo(valormax,i,a,b)[0], MonteCarlo(valormax,i,a,b)[1]
        valores_I.append(I)
        desv_valores.append(desv)
 

plt.title('Valores de la integral respecto al número de pasos y su desviación') 
plt.errorbar(n_pasos,valores_I, yerr = desv_valores, fmt = '-o', ecolor = 'orange')
plt.ylabel('Valores de la integral')
plt.xlabel('Nº de pasos')
plt.show()    

#4

intervalo_2 = np.linspace(0,2,len(x_dentro))
plt.plot(x,f(x))
plt.title('Representación de los puntos obtenidos con el método MonteCarlo')
plt.plot(x_dentro,y_dentro,'o',label = 'Valores dentro')
plt.plot(x_fuera,y_fuera,'o', label = 'Valores fuera')
plt.ylim(0,1.)

plt.legend(loc = 'best')
plt.show()

#EXTRA  1
#función error

def f_2(x):
    return (np.exp(-x**2))


x = np.linspace(-3,3,nt)
a = x[0]
b = x[-1]

plt.plot(x, f_2(x))
plt.title('Valor de la función entre $x = -3$ y $x = 3$')
plt.xlabel('Número de puntos')
plt.ylabel('f(x)')
plt.show()
plt.close()
#probablemente es posible generar una única función MonteCarlo que sea capaz de aceptar cualquier f(x)
#pero no tengo claro si siguiendo esta lógica es posible
def MonteCarlo2(valormax,nt,a,b):
    y = np.random.uniform(0,valormax,nt)
    x = np.random.uniform(a,b,nt) 
 
    x_dentro = []
    y_dentro = []
    x_fuera = []
    y_fuera = []
    for i in range(0,nt):
        valor = f_2(x[i])
        if y[i] < valor:
            x_dentro.append(x[i])
            y_dentro.append(y[i])
        else:
            x_fuera.append(x[i])
            y_fuera.append(y[i])
   
    x_dentro,y_dentro = np.array(x_dentro),np.array(y_dentro)
    frac_pares = len(x_dentro)/nt
    I = valormax*(b-a)*frac_pares       
    desv = np.sqrt(I*(valormax*(b-a)-I))/np.sqrt(nt)
    return I,desv,x_dentro,y_dentro,x_fuera,y_fuera

resultados = MonteCarlo2(valormax,nt,a,b)
I,desv,x_dentro,y_dentro,x_fuera,y_fuera = resultados
print(f'El valor de la integral es {I:.3f} ± {desv:.3f}')

valores_I = []
n_pasos = []
desv_valores = []

if nt <= 5000:
    for i in range(1,nt,10):
        n_pasos.append(i)
        I,desv = MonteCarlo(valormax,i,a,b)[0], MonteCarlo(valormax,i,a,b)[1]
        valores_I.append(I)
        desv_valores.append(desv)
if nt >= 100000:
    for i in range(1,nt,10000):
        n_pasos.append(i)
        I,desv = MonteCarlo(valormax,i,a,b)[0], MonteCarlo(valormax,i,a,b)[1]
        valores_I.append(I)
        desv_valores.append(desv)
else:
    for i in range(1,nt,1000):
        n_pasos.append(i)
        I,desv = MonteCarlo(valormax,i,a,b)[0], MonteCarlo(valormax,i,a,b)[1]
        valores_I.append(I)
        desv_valores.append(desv)
 

plt.title('Valores de la integral respecto al número de pasos y su desviación') 
plt.errorbar(n_pasos,valores_I, yerr = desv_valores, fmt = '-o', ecolor = 'orange')
plt.ylabel('Valores de la integral')
plt.xlabel('Nº de pasos')
plt.show()    

plt.plot(x,f_2(x))
plt.title('Representación de los puntos obtenidos con el método MonteCarlo')
plt.plot(x_dentro,y_dentro,'o',label = 'Valores dentro')
plt.plot(x_fuera,y_fuera,'o', label = 'Valores fuera')
plt.ylim(0,1.)

plt.legend(loc = 'best')
plt.show()

#EXTRA 2
x = np.linspace(-1,1,nt)
y = np.linspace(-1,1,nt)
a = x[0]
b = x[-1]
valormax = 1
def f_3(x,y):
    return (x)**2+y**2 
    
def MonteCarloCirculo(valormax,nt,a,b):
    y = np.random.uniform(-valormax,valormax,nt)
    x = np.random.uniform(a,b,nt) 
    
    x_dentro = []
    y_dentro = []
    x_fuera = []
    y_fuera = []
    for i in range(0,nt):
        valor = f_3(x[i],y[i])
        if 1 > valor:
            x_dentro.append(x[i])
            y_dentro.append(y[i])
        else:
            x_fuera.append(x[i])
            y_fuera.append(y[i])
    x_dentro,y_dentro = np.array(x_dentro),np.array(y_dentro)
    frac_pares = len(x_dentro)/nt
    A = (b-a)*2*valormax*frac_pares
    return A,x_dentro,y_dentro,x_fuera,y_fuera

resultados = MonteCarloCirculo(valormax,nt,a,b)
A,x_dentro,y_dentro,x_fuera,y_fuera = resultados

print(f'Valor del área es {A:.3f}')
#plotamos el círculo a ver qué sale
plt.title('Representación de los puntos obtenidos con el método MonteCarlo')
plt.plot(x_dentro,y_dentro,'o', label = 'puntos dentro', color = 'orange')
plt.plot(x_fuera,y_fuera,'o',label = 'puntos fuera',color = 'green')
plt.legend(loc = 'best')
plt.show()
#creamos un array para calcular el área en todas las dimensiones
N = 10
dim = []

a = -1
b = 1
def MonteCarloHiperesfera(nt,a,b):
    puntos = np.random.uniform(-1,1,(nt,10))
    puntos_dentro = []
    puntos_fuera = []
    for i in range(nt):
        norma = np.linalg.norm(puntos[i])
        if 1 >= norma:
            puntos_dentro.append(puntos[i])
        else:
            puntos_fuera.append(puntos[i])
    puntos_dentro, puntos_fuera = np.array(puntos_dentro),np.array(puntos_fuera)
    frac_pares = len(puntos_dentro)/nt
    V = 2**10*frac_pares
    return V
    
total_vol = MonteCarloHiperesfera(nt,a,b)
print(f'El volumen total de la hiperesfera en {N} dimensiones es {total_vol:.3f}')    
