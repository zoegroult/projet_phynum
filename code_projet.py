import numpy as np 
import matplotlib.pyplot as plt
import code_preprojet as cpp


# variables globales
H0 = 0.070      # 67.4*1e-19*1e9*365.25*24*3600/3.03 : cste de Hubble à notre ère (conversion en Gyr-1)
ai = 1e-6     # facteur d'échelle initial (au début de l'Univers) 
ti = 1/(365.25*1e9)   # i = temps initial (vers les débuts de l'Univers) en Gyr
t0 = 14 # 0 = temps de notre ère (age de l'Univers) en Gyr
G = 6.6742*1e-11  # constante gravitationnelle de Newton ( À convertir ???)


Omega_m0 = 0.315
Omega_r0 = 9*1e-5
Omega_L0 = 0.685


mu_x = mu_y = 0
sigma_x = sigma_y = 1
rho = 

A,T = cpp.resolution_EDO_Friedmann()

Omega_L = cpp.Omega_l(A)
Omega_M = cpp.Omega_m(A)
Omega_R = cpp.Omega_r(A)

rho_C = cpp.rho_c(A)

rho_L = Omega_L*rho_C
rho_R = Omega_R*rho_C
rho_M = Omega_M*rho_C


#------------------------------------------------------------------------------------------------------------------------------------
# partie servant a trouver le taux de croissance linéaire D(a)
def f(x):
    return (Omega_m0 * x**-1  + Omega_L0 * x**2  + Omega_r0)**(-3/2)


def H(x):
    return H0 * np.sqrt(Omega_r0 * x**-4  +  Omega_m0 * x**-3  +  Omega_L0)


def integration_rectangle(a, N=1e4):
    Somme = np.zeros( len(a) )

    for i in range( len(a) ):
        print(i)
        x = np.linspace(1e-7, a[i], int(N))
        
        S = 0
        S += 0.5 * (x[1:]-x[:-1]) * (f(x[1:]) + f(x[:-1]))
        
        #for j in range( len(x)-1 ):
        #    S += 0.5 * (x[j+1] - x[j]) * (f(x[j+1]) + f(x[j]))

        Somme[i] = S  
    return Somme   


def D(a):
    return integration_rectangle(a) * H(a)/H0 






#------------------------------------------------------------------------------------------------------------------------------------
# Partie servant a trouver le champ de potentiel ∇Ψ :

def phi_i(x, y, sigma_x, sigma_y, mu_x, mu_y, rho ):    #rho est la corrélation entre x et y 
    X = x - mu_x
    Y = y - mu_y
    SIGMA = np.array([[sigma_x**2, rho*sigma_x*sigma_y] , [rho*sigma_x*sigma_y, sigma_y**2]])

    return (2*np.pi*np.sqrt(np.linalg.det(SIGMA)))**-1 * np.exp(-0.5 * (np.linalg.inv(SIGMA)[0,0] * X**2  + X*Y*(np.linalg.inv(SIGMA)[0,1] + np.linalg.inv(SIGMA)[1,0]) +  np.linalg.inv(SIGMA)[1,1] * Y**2))
    


def psi(x, y, sigma_x, sigma_y, mu_x, mu_y, rho ):
    return -2*phi_i(x, y, sigma_x, sigma_y, mu_x, mu_y, rho ) / (3*H0*Omega_m0)




def nabla_psi(x,y, sigma_x, sigma_y, mu_x, mu_y, rho):
    
    PSI = psi(x,y, sigma_x, sigma_y, mu_x, mu_y, rho)

    dx = np.abs(X[0,0]-X[0,1])
    dy = np.abs(Y[0,0]-Y[1,0])

    dpsi_dx[1:-1, :] = (PSI[2:,:] - PSI[:-2, :]) / (2*dx)
    dpsi_dy[:, 1:-1] = (PSI[:,2:] - PSI[:, :-2]) / (2*dy)

    # les bords :
    dpsi_dx[0,:] = (PSI[1,:] - PSI[0,:])/ dx
    dpsi_dx[-1,:] = (PSI[-1,:] - PSI[-2,:])/ dx

    dpsi_dy[:,0] = (PSI[:,1] - PSI[:,0])/ dy
    dpsi_dx[:,-1] = (PSI[:,-1] - PSI[:,-2])/ dy

    return dpsi_dx, dpsi_dy
    




#------------------------------------------------------------------------------------------------------------------------------------
# Partie servant à avoir la position x de chaque particules en fonction du temps :

def position(x, y, a, t):
    X, Y = np.meshgrid(x,y)
    dpsi_dx, dpsi_dy = nabla_psi(X, Y, sigma_x, sigma_y, mu_x, mu_y, rho)

    xx = X[0,:] + D(a) * dpsi_dx[0,:]
    yy = 

    



plt.figure(figsize=(7,5))
plt.plot( cpp.A, D, 'bo:')
plt.title("évolution du taux de croissance linéaire")
plt.xlabel("facteur d'échelle a(t)")
plt.ylabel("taux de croissance linéaire D(a)")
plt.show()
