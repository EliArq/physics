# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#Randomwalk method

import numpy as np
import matplotlib.pyplot as plt

N = int(input('Número de pasos: '))
I = int(input('Número de iteraciones: '))

Pr = float(np.random.random_sample())

Pl = 1-Pr
#print(Pr,Pl)
longitud = np.random.rand(2)

Nr = np.random.randint(0,N+1)

Nl = N - Nr

def Random_walk(N, Pr,Pl,longitud):
 
    x = np.zeros(N+1)

    for i in range(1,N+1):
        salto = np.random.choice(['L','R'],p=[Pl,Pr])
        
        if salto == 'L':
            x[i] = x[i-1] - longitud[0]
            
        else:
            x[i] = x[i-1] + longitud[1]
    return x



puntos = []

for i in range(I):
    x = (Random_walk(N, Pr,Pl,longitud))
    puntos.append(x[-1])


plt.title('Random walk $n = '+str(N) + '$ saltos,' + str(I) + ' interacciones')    
plt.hist(puntos, bins = 30, density = True, histtype='barstacked')
plt.xlabel('Pasos')
plt.ylabel('P(x)')
plt.show()
plt.close()

media = np.mean(puntos)
desviacion = np.std(puntos)
norma_L = np.sqrt(longitud[0]**2+longitud[1]**2)
prob_Nr = N*Pr
prob_Nl = N*Pl


media_teo = N*(Pr*longitud[1] - Pl*longitud[0])

media_teo2 = ((longitud[0]+longitud[1])**2)*((Pr**2)*(N**2)+N*Pr*Pl)-2*N*((longitud[0]**2)+longitud[1]*longitud[0])*(N*Pr)+(longitud[0]**2)*(N**2)
var_teo = media_teo2-media_teo**2
desv_teo = np.sqrt(var_teo)
print('Media teórica: ', media_teo)
print('Desviación teórica: ', desv_teo)
print('La media de la distribución obtenida con numpy es: ', media, 'mientras que la desviación es: ', desviacion)

def gauss(x, media,desviacion):
    
    return 1./(np.sqrt(np.pi*2.)*desviacion)*np.exp(-(x-media)**2/(2*desviacion**2))

x_gauss = np.linspace(min(puntos),max(puntos),I)

plt.title(' Distribución gaussiana')
plt.plot(x_gauss, gauss(x_gauss,media,desviacion))
plt.stem(x_gauss, gauss(x_gauss,media,desviacion), linefmt='blue', markerfmt='go', basefmt=" ", label = 'Valores de la distribución')
plt.xlabel('Pasos')
plt.ylabel('P(x)')
plt.legend()
plt.show()
plt.close()

lambda_pois = N*Pr
m = Nr-Nl
eventos = np.copy(puntos)


def Poisson(Nr,lambda_pois):
    try:
        return lambda_pois**Nr*np.exp(-lambda_pois)/(np.math.factorial(Nr))
    except:
         return 0

poisson_vals = []
prob_poisson = np.zeros(len(eventos))

for Nr in range(0,len(eventos)):

    poisson_vals.append(Poisson(Nr,lambda_pois))
    prob_poisson[Nr] = np.array(Poisson(Nr,lambda_pois))

prob_norm =  prob_poisson/sum(prob_poisson)


xtotal = Random_walk(N,Pr,Pl,longitud)

poisson_muestra = np.random.choice(eventos, len(xtotal),True, p= prob_norm)

plt.title('Distribución Poisson')
plt.plot(x_gauss,poisson_vals)
plt.stem(x_gauss, poisson_vals, linefmt='green', markerfmt='go', basefmt=" ", label = 'Valores de la distribución')
plt.xlabel('Pasos')
plt.ylabel('P(x)')
plt.legend()
plt.show()

plt.title('Random walk $n = '+str(N) + '$ saltos,' + str(I) + ' interacciones') 
#plt.plot(x_gauss, gauss(x_gauss,media,desviacion), label = 'Gauss') 
plt.hist(puntos, bins = 30, density = True, histtype='barstacked', color = 'green', label = 'Randomwalk')
plt.hist(poisson_muestra, bins = 30, density = True, color = 'orange', label = 'Poisson')
plt.xlabel('Pasos')
plt.ylabel('P(x)')
plt.legend()

plt.show()
plt.close()

