#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 17:00:15 2025

@author: zoegroult
"""

import numpy as np 
import matplotlib.pyplot as plt
import math 


# variables globales
H0 = 0.070                          # 67.4*1e-19*1e9*365.25*24*3600/3.03 : cste de Hubble à notre ère (conversion en Gyr-1)
ai = 1e-6                           # facteur d'échelle initial (au début de l'Univers) 
ti = 1/(365.25*1e9)                 # i = temps initial (vers les débuts de l'Univers) en Gyr
t0 =  14                            # 0 = temps de notre ère (age de l'Univers) en Gyr
G = 6.6742*1e-11                    # constante gravitationnelle de Newton ( À convertir ???)
Omega_m0 = 0.315
Omega_r0 = 9*1e-5
Omega_L0 = 0.685



def resolution_complète_EDO_Friedmann(om_0=Omega_m0, or_0=Omega_r0, ol_0=Omega_L0, h=1e-5):   # fonction qui resoud l'equa diff avec la methode d'euler (h=pas de temps) 
    N = int((t0-ti)/h)
    
    A=np.zeros(N+1)
    A[0]=ai
    
    T=np.zeros(N+1)
    T[0]=ti
    for i in range(N):
        A[i+1]=A[i] + h * H0 * np.sqrt( (A[i]**-2) * or_0 + (A[i]**-1) * om_0 + (A[i]**2) * ol_0 )
        T[i+1]=T[i] + h
        
    return A, T




def a(t_voulu, om_0=Omega_m0, or_0=Omega_r0, ol_0=Omega_L0, H = 1e-5):
    """kljcbsqmkvbzqmjbvmqjbvkbvmkqjzbvmozhbvqjbzrlvbzmjoebgvmouzbvzmjbv"""

    print('a')
    A = ai
    t = ti

    while t < t_voulu :

        if A < 1e-5:
            h = 1e-12

        elif A < 1e-2:
            h = 1e-10
   
        else:
            h = H

        A = A + h * H0 * np.sqrt( (A**-2) * or_0 + (A**-1) * om_0 + (A**2) * ol_0 )
        t += h

    return A





def a2(t_voulu, om_0=Omega_m0, or_0=Omega_r0, ol_0=Omega_L0, H = 1e-5, afficher = False):
    """kljcbsqmkvbzqmjbvmqjbvkbvmkqjz
    
    bvmozhbvqjbzrlvbzmjoebgvmouzbvzmjbv"""
    A = ai
    t = ti

    while t < t_voulu :

        if t < 1e-10:
            h = 1e-12

        elif t < 1e-8:
            h = 1e-10

        elif t < 1e-6:
            h = 1e-8

        elif t< 1e-4:
            h = 1e-6

        else:
            h = H


        A = A + h * H0 * np.sqrt( (A**-2) * or_0 + (A**-1) * om_0 + (A**2) * ol_0 )
        t += h
    
    if afficher == True:
        print(f'a({t_voulu}) =', A/1.0154729584922937)
        
    return A/1.0154729584922937







def Omega_m(t, om_0=Omega_m0, or_0=Omega_r0, ol_0=Omega_L0):  # parametre de densité de la matiere totale (mat visible + mat noire)
    
    Omega_m = om_0 / ( (a(t)**-1)*or_0 + om_0 + (a(t)**3) * ol_0 )
    
    return Omega_m





def Omega_r(t, om_0=Omega_m0, or_0=Omega_r0, ol_0=Omega_L0):  # parametre de densité de radiation

    Omega_r = or_0 / ( or_0 + a(t) * om_0 + (a(t)**4) * ol_0 ) 
    
    return Omega_r





def Omega_l(t, om_0=Omega_m0, or_0=Omega_r0, ol_0=Omega_L0):     # parametre de densité d'énergie noire 

    Omega_l = ol_0 / ( (a(t)**-4) * or_0 + (a(t)**-3) * om_0 + ol_0 ) 
    
    return Omega_l





def rho_c(t, om_0=Omega_m0, or_0=Omega_r0, ol_0=Omega_L0):

    rho_c = (3/(8*np.pi*G)) * H0**2 * ((a(t)**-4) * or_0 + (a(t)**-3) * om_0 + ol_0)
    
    return  rho_c




"""lines = []
with open("A_T_Om_Ol_Or_rhoc.dat", 'w') as f:
    for j in range(len(A)):
        lines.append(f"{A[j]} {T[j]} {Omega_M[j]} {Omega_L[j]} {Omega_R[j]} {rho_critique[j]}\n")
    f.writelines(lines)

data = np.loadtxt("A_T_Om_Ol_Or_rhoc.dat")
A = data[:,0]
T = data[:,1]
Omega_L = data[:,2]
Omega_M = data[:,3]
Omega_R = data[:,4]
rho_C = data[:,5]



# plot de a(t)
t = np.linspace(1e-7, 14, 1000)
A = []
for i in t:
    A.append(a2(i))
plt.figure(figsize=(7,5))
plt.plot(t, A, 'b:')
plt.xlabel('temps t en Gyr')
plt.ylabel('facteur d\'échelle a(t)')
plt.title(f"Courbe d\'évolution du facteur d\'échelle a(t),  t0={t0}")
plt.grid(True)
plt.show()




# plot des Omega_i(t)
plt.figure(figsize=(7,5))
plt.plot(T, Omega_L, 'r:', label="Omega_L(t)")
plt.plot(T, Omega_M, 'b:', label="Omega_M(t)")
plt.plot(T, Omega_R, 'g:', label="Omega_R(t)")
#plt.xscale("log")
plt.xlabel("temps t (en Gyr)")
plt.ylabel("Omega_i(t)")
plt.title(f"Évolution des paramètres de densité en fonction du temps , t0={t0}")
plt.legend()
plt.grid(True)
plt.show()



# plot des Omega_i(a) en log/log: 
plt.figure(figsize=(7,5))
plt.plot(A, Omega_L, 'r:', label="Omega_L(t)")
plt.plot(A, Omega_M, 'b:', label="Omega_M(t)")
plt.plot(A, Omega_R, 'g:', label="Omega_R(t)")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("facteur d'échelle a(t)")
plt.ylabel("Omega_i(t)")
plt.title(f"Évolution des paramètres de densité en fonction du facteur d'échelle, t0={t0}")
plt.legend()
plt.grid(True)
plt.show()



# plot des rho_i(a) en log/log:
plt.figure(figsize=(7,5))
plt.plot(A, Omega_L*rho_C, 'r:', label="rho_L(t)")
plt.plot(A, Omega_M*rho_C, 'b:', label="rho_M(t)")
plt.plot(A, Omega_R*rho_C, 'g:', label="rho_R(t)")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("facteur d'échelle a(t)")
plt.ylabel("rho_i(t)")
plt.title(f"Évolution des densités en fonction du facteur d'échelle en log/log, t0={t0}")
plt.legend()
plt.grid(True)
plt.show()




# plot des rho_i(t) en log/log:
plt.figure(figsize=(7,5))
plt.plot(T, Omega_L*rho_C, 'r:', label="rho_L(t)")
plt.plot(T, Omega_M*rho_C, 'b:', label="rho_M(t)")
plt.plot(T, Omega_R*rho_C, 'g:', label="rho_R(t)")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("temps t en Gyr")
plt.ylabel("rho_i(t)")
plt.title(f"Évolution des densités en fonction du temps en log/log, t0={t0}")
plt.legend()
plt.grid(True)
plt.show()"""



