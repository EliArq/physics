# -*- coding: utf-8 -*-
"""
Created on Sat Dec  7 11:38:08 2024

@author: Usr
"""
#Fourier transform using different methods

import scipy.fft as scifft 
import numpy as np
import matplotlib.pyplot as plt
import pygame
import os
import librosa
import time


y, sr = librosa.load(librosa.example('trumpet'))

print(sr)
tempo, beats = librosa.beat.beat_track(y = y, sr = sr)

t = librosa.frames_to_time(beats,sr=sr) #np.linspace de tiempos del audio
print(t)
print(tempo)

#DFT


start = time.time()
def DFT(x):
    N = len(x)
    n = np.arange(N)
    k = n.reshape((N,1))
    e = np.exp(-2j*np.pi*k*n/N)
    X = e@x
    
    return X


N = len(y)
fragment_t = 1024
n_frag = len(y)//fragment_t

print(n_frag)

dftprom = []
frecuencias = []

#en esta parte se va a coger solo las frecuencias positivas, porque por simetría tienes duplicación de datos
#he quitado que las frecuencias se obtengan con numpy.fft, no nos interesa para hacerlo “a mano”, se usa una fórmula que te relaciona los “batidos” para cada fragmento de onda y se divide entre el número total

for i in range(n_frag):
    fragmento = y[i*fragment_t:(i+1)*fragment_t]
    if len(fragmento) == 1:
        mitad_N = n_frag//2
    else:
        mitad_N = len(fragmento)//2
    dft = DFT(fragmento)/len(fragmento)
    amplitudes = np.abs(dft[:mitad_N])
    dftprom.append(amplitudes)

dft_tiempo = time.time() - start
print(f'Tiempo de ejecución DFT: {dft_tiempo:.3f}')    


dftprom = np.array(dftprom)
if len(dftprom) > 0:
    promedio_dft = np.mean(dftprom,axis = 0) 
else:
    promedio_dft = dftprom[0]
freq = np.linspace(0,sr/2,len(promedio_dft),endpoint = False)
#frecuencias = np.linspace(0,sr/2,len(dftprom),endpoint = False)
freq_dominante_dft = freq[np.argmax(promedio_dft)]
print(f'Frecuencia dominante en DFT: {freq_dominante_dft} Hz')
print('Periodos de señal', len(freq)/2)
plt.figure(figsize =(12,4))
tiempo = np.linspace(0,len(y)/sr,len(y))
plt.plot(tiempo, y)
plt.title('Señal en el dominio de tiempo')
plt.xlabel('Tiempo (s)')
plt.ylabel('Amplitud')

plt.show()

#esto es lo más interesante, la gráfica de la transformada
#si usamos la partición de frecuencias merece más la pena hacer el promedio, pero con la trompeta coge todo el espectro

plt.figure(figsize=(12, 4))
plt.stem(freq, promedio_dft, markerfmt=" ", basefmt="-b")
plt.title("Espectro de frecuencias (DFT)")
plt.xlabel("Frecuencia (Hz)")
plt.ylabel("DFT Amplitud |Y(frecuencia)|")
plt.show()

#FFT
#voy a probar con el mismo número de fragmentos

frecuencias2 = []
start = time.time()
def FFT(x):
    N = len(x)
    if N == 1:
        return x
    else:
        X_even = FFT(x[::2])
        X_odd = FFT(x[1::2])
        factor = np.exp(-2j*np.pi*np.arange(N)/N)
        X = np.concatenate([X_even+factor[:int(N/2)]*X_odd,X_even+factor[int(N/2):]*X_odd])
        return X


fftprom = []
for i in range(n_frag):
    fragmento = y[i*fragment_t:(i+1)*fragment_t]

    if len(fragmento) == 1:
        mitad_N = n_frag//2
    else:
        mitad_N = len(fragmento)//2
    fft = FFT(fragmento)/len(fragmento)
    amplitudes = np.abs(fft[:mitad_N])
    fftprom.append(amplitudes)
 
fft_tiempo = time.time() - start
print(f'Tiempo de ejecución FFT: {fft_tiempo:.3f}') 


#fft = FFT(y)
fftprom = np.array(fftprom)
if len(fftprom) > 0:
    promedio_fft = np.mean(fftprom,axis = 0) 
else:
    promedio_fft = fftprom[0]
    
freq = np.linspace(0,sr/2,len(promedio_fft),endpoint = False)

freq_dominante_fft = freq[np.argmax(promedio_fft)]
print(f'Frecuencia dominante en FFT: {freq_dominante_fft} Hz')

#frecuencias2 = np.linspace(0,sr/2,len(dftprom),endpoint = False)
plt.figure(figsize = (12,4))
plt.stem(freq,promedio_fft,markerfmt = ' ',basefmt = "-b")
plt.title('Espectro de frecuencias FFT')
plt.xlabel('Frecuencia (Hz)')
plt.ylabel('FFT Amplitud |Y(freccuencia)|')
plt.show()

#IDFT

start = time.time()
def IDFT(X):
    N = len(X)
    n = np.arange(N)
    k = n.reshape((N,1))
    e = np.exp(2j*np.pi*k*n/N)
    x = (1/N)*(e@X)

    return x


idftprom = []
for i in range(n_frag):
    fragmento = y[i*fragment_t:(i+1)*fragment_t]
    if len(fragmento) == 1:
        mitad_N = n_frag//2
    else:
        mitad_N = len(fragmento)//2
    dft = DFT(fragmento)
    idft = IDFT(dft).real
    idftprom.append(idft)

reconstruida = np.concatenate(idftprom)
idft_tiempo = time.time() - start
print(f'Tiempo de ejecución IDFT: {idft_tiempo:.3f}') 
idftprom = np.array(idftprom)
if len(idftprom) > 0:
    promedio_idft = np.mean(idftprom,axis = 0) 
else:
    promedio_idft = idftprom[0]
freq = np.linspace(0,sr/2,len(promedio_idft),endpoint = False)

#frecuencias3 = np.linspace(0,sr/2,len(dftprom),endpoint = False)
tiempo = np.linspace(0,len(reconstruida)/sr, len(reconstruida))
plt.figure(figsize = (12,4))
plt.title('Señal reconstruida IDFT')
plt.plot(tiempo,reconstruida)
plt.xlabel('Tiempo (s)')
plt.ylabel('IDFT Amplitud')
plt.show()

#USANDO SCIPY

#FFT
start = time.time()
FFTscipy = scifft.fft(y, axis = 0)
fftsci_tiempo = time.time() - start
print(f'Tiempo de ejecución FFT Scipy: {fftsci_tiempo:.3f}') 
ts = 1/(sr)  #pasos de tiempo
print(ts)
freq = scifft.fftfreq(N, d =ts)
amplitudes = np.abs(FFTscipy)
freq_dominante_scipy = freq[np.argmax(amplitudes)]
print(f'Frecuencia dominante en FFT de Scipy: {freq_dominante_scipy} Hz')
print('Periodos de señal', len(freq)/2)

plt.figure(figsize = (12,4))
plt.stem(freq[1:N//2],np.abs(FFTscipy[1:N//2])/len(y), markerfmt = ' ',basefmt = "-b")
plt.title('FFT Scipy')
plt.xlabel('Frecuencia (Hz)')
plt.ylabel('FFT Amplitud |Y(freccuencia)|')
plt.show()

#IDFT

start = time.time()
IFFTscipy = scifft.ifft(FFTscipy).real
ifft_tiempo = time.time() - start
print(f'Tiempo de ejecución IFFT: {ifft_tiempo:.3f}') 
tiempo2 = np.linspace(0,len(y)/sr, len(IFFTscipy))
plt.figure(figsize = (12,4))
plt.plot(tiempo2, IFFTscipy)
plt.title('Señal reconstruida IFFT Scipy')
plt.xlabel('Tiempo (s)')
plt.ylabel('IFFT Amplitud ')
plt.show()


