#Nuclear phyics script to obtain Semi-empirical mass formula and binding energy in Uranium-238 and U-235 besides isotopes
import matplotlib.pyplot as plt
import numpy as np
import os

#Simple comprobation for changing path of working directory
#current_dir = os.getcwd()
#change = os.chdir(os.path.join(current_dir," "))
#print("Current Directory:", os.getcwd())
#Semi-empirical masss formula
def B(A,Z):
    if A%2 == 0 and Z%2 ==0:
        print('Número A: ',A,' y Z: ',Z,' término paridad: ',11.2/np.sqrt(A))
        return (15.8*A - 18.3*A**(2/3) - 0.714*(Z**2)/(A**(1/3)) - 23.2*(A - 2*Z)**2/A + 11.2/np.sqrt(A))
        
    else:
        if A%2 != 0 and Z%2 != 0:
            print('Número A: ',A,' y Z: ',Z,' término paridad: ',-11.2/np.sqrt(A))
            return (15.8*A - 18.3*A**(2/3) - 0.714*(Z**2)/(A**(1/3)) - 23.2*(A - 2*Z)**2/A - 11.2/np.sqrt(A))
        else:
            return (15.8*A - 18.3*A**(2/3) - 0.714*(Z**2)/(A**(1/3)) - 23.2*(A - 2*Z)**2/A)

Z = 92
A = [235,238]

for i in A:
    Energia = B(i,Z)/i
    if Z**2/i >= 50:
        print('Cumple la condición de fisión', Z**2/i)
    print(f"Energia para el Uranio {i:.2f} : ", Energia, 'MeV')


#Isotopes
A2 = np.linspace(214,240,36)

Energia2 = []
for i in A2:
    Energia2.append(B(i,Z)/i)


plt.plot(A2,Energia2, '.-')
plt.xlabel('A')
plt.ylabel('Energía por nucleón (MeV)')
plt.title('Energía de enlace por nucleón para isótopos del Uranio')
#plt.savefig('Energia_Uranio.png')
plt.show()

#Energy values 

Particulas = {
    'Uranio-238': (238,92),
    'Kriptón-92': (92,36),
    'Bario-141': (141,56)
}

Energia3 = []
for i in Particulas:
    nombres = i[0:10]
    A,Z = Particulas[i]
    Energia3.append(B(A,Z))
    print(f'Energia para {nombres}:', Energia3[-1], 'MeV')


A = [Particulas[i][0] for i in Particulas]
Z = [Particulas[i][1] for i in Particulas]
Q = -Energia3[0] + (Energia3[1] + Energia3[2])
Qprods = Energia3[1] + Energia3[2]
print('Energía liberada en la fisión:', Q, 'MeV')

nucleos = ['U-238', 'Kr-92', 'Ba-141', 'Productos\nCombinados', 'Energía\nTotal']
energias = [Energia3[0], Energia3[1], Energia3[2], Qprods, Q]

plt.bar(nucleos, energias, color=['blue', 'orange', 'green', 'red', 'purple'])
plt.ylabel('Energía de enlace (MeV)')
plt.title('Energía de enlace para productos de fisión del Uranio')
#plt.savefig('Energia_Fision_Uranio.png')
plt.show()

