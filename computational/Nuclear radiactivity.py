#Script of nuclear physics to obtain radiactivity of U-235 and U-238
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import os


#current_dir = os.getcwd()
#change = os.chdir(os.path.join(current_dir,""))
#print("Current Directory:", os.getcwd())

t = np.linspace(0, 10**9, 1000)
#Number and activity of particles over time
def N(t):
    Tmedia = [4.468*10**9, 703.8*10**6] #U-238, U-235 en años
    lmbd = []
    for i in Tmedia:
        lmbd.append(np.log(2)/i)
    N0 = 1*10**6
    N = np.zeros((2, len(t)))
    A = np.zeros((2, len(t)))
    indices = np.linspace(0, len(t)-1, 15, dtype=int)  
    for i in range(len(t)):
        N[:,i] = N0*np.exp(-(np.dot(lmbd, t[i])))
        A[:,i] = lmbd*N[:,i]
        if i in indices:
            print(f'Número de partículas U-238 y U-235 después de {t[i]} años: {N[:,i]}')
            print(f'Actividad del núcleo después de {t[i]} años: {A[:,i]} 1/años')
    return N, A, Tmedia, lmbd


N, A, Tmedia, lmbd = N(t)

plt.plot(t,N[0,:], label = 'Núcleos de U-238')
plt.plot(t,N[1,:], label = 'Núcleos de U-235')
plt.legend()
plt.yscale('log')
plt.title('Decaimiento de los núcleos radiactivos')
plt.xlabel('Tiempo (Años)')
plt.ylabel('Número de núcleos restantes')
plt.tight_layout() 
#plt.savefig('Decaimiento_Uranio.png')
plt.show()


plt.plot(t,A[0,:], label = 'Actividad de U-238')
plt.plot(t,A[1,:], label = 'Actividad de U-235')
plt.legend()
plt.yscale('log')
plt.title('Actividad de los núcleos radiactivos')
plt.xlabel('Tiempo (Años)')
plt.ylabel('Actividad (Desintegraciones por año)')
plt.tight_layout() 
#plt.savefig('Actividad_Uranio.png')
plt.show()

#First decay chain
Tmedia.append(24.1/365) #Th234
Tmedia.append(25.52/(24*365)) #Th231
def decaimiento(Tmedia):
    lmbd = np.zeros(len(Tmedia))
    for i in range(len(Tmedia)):
        lmbd[i] = np.log(2)/Tmedia[i]

    def sistema_ecuaciones(t,N2):
        dNdt = np.zeros(4) 
        #Cadena u-238 -> th-234
        dNdt[0] = -lmbd[0]*N2[0]
        dNdt[2] = -lmbd[2]*N2[2] + lmbd[0]*N2[0]
        #Cadena u-235 -> th-231
        dNdt[1] = -lmbd[1]*N2[1] 
        dNdt[3] = -lmbd[3]*N2[3] + lmbd[1]*N2[1]

        return dNdt

    cond_iniciales = [1*10**6, 1*10**6, 0, 0]
    sol_ec = solve_ivp(sistema_ecuaciones, (t[0], t[-1]), cond_iniciales, t_eval = t, method = 'BDF').y
    N2 = sol_ec
    A2 = np.zeros((4, len(t)))
    for i in range(len(N2)):
        A2[i,:] = lmbd[i]*N2[i,:]

    return N2, A2,lmbd


N2, A2,lmbd = decaimiento(Tmedia)
indices = np.linspace(0, len(t)-1, 15, dtype=int)  
for i in range(len(t)):
    if i in indices:
        print(f'Número de partículas + hijas después de {t[i]} años: {N2[:,i]}')
        print(f'Actividad de todas las particulas después de {t[i]} años: {A2[:,i]} 1/años')

print(lmbd)

plt.plot(t,N2[0,:], label = 'Núcleos de U-238')
plt.plot(t,N2[1,:], label = 'Núcleos de U-235')
plt.plot(t,N2[2,:], label = 'Núcleos de Th-234')
plt.plot(t,N2[3,:], label = 'Núcleos de Th-231')
plt.legend()
plt.yscale('log')
plt.title('Decaimiento de los núcleos radiactivos')
plt.xlabel('Tiempo (Años)')
plt.ylabel('Número de núcleos restantes')
plt.tight_layout() 
#plt.savefig('Decaimiento_Uranio_hijas.png')
plt.show()

plt.plot(t,A2[0,:], label = 'Actividad de U-238')
plt.plot(t,A2[1,:], label = 'Actividad de U-235')
plt.plot(t,A2[2,:], label = 'Actividad de Th-234')
plt.plot(t,A2[3,:], label = 'Actividad de Th-231')
plt.legend()
plt.yscale('log')
plt.title('Actividad de los núcleos radiactivos')
plt.xlabel('Tiempo (Años)')
plt.ylabel('Actividad (Desintegraciones por año)')
plt.tight_layout() 
#plt.savefig('Actividad_Uranio_hijas.png')
plt.show()


# Calculations for 90%, 50% and 10% of mass decay over time
N0 = 1*10**6
decae0_9 = 0.9*N0
decae0_5 = 0.5*N0
decae0_1 = 0.1*N0

particulas = {
    'U-238': N2[0,:],
    'U-235': N2[1,:],
    'Th-234': N2[2,:],
    'Th-231': N2[3,:],

}

for k, nucleo in enumerate(particulas.keys()):
    print(f"\n{nucleo}:")
    vida_media = Tmedia[k]
    lambda_k = lmbd[k]
    
    for porcentaje, etiqueta in zip([90, 50, 10], ['10%', '50%', '90%']):
        
        tiempo_teorico = -np.log(porcentaje/100) / lambda_k
        print(f"  Decaimiento del {etiqueta}: {tiempo_teorico:.2e} años")
        print(f"    Relación con vida media: {tiempo_teorico/vida_media:.2f}")
        
        if porcentaje == 50:
            print(f"    Vida media teórica: {vida_media:.2e} años")