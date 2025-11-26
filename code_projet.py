import numpy as np 
import matplotlib.pyplot as plt
import code_preprojet as cpp
from math import *



# variables globales
H0 = 0.070      # 67.4*1e-19*1e9*365.25*24*3600/3.03 : cste de Hubble à notre ère (conversion en Gyr-1)
ai = 1e-6     # facteur d'échelle initial (au début de l'Univers) 
ti = 1/(365.25*1e9)   # i = temps initial (vers les débuts de l'Univers) en Gyr
t0 = 14 # 0 = temps de notre ère (age de l'Univers) en Gyr

M_soleil = 1.98847e30  # masse du soleil en kg
Gyr_en_s = 3.15576e16  # 1 Gyr = 3.15576e16 s
Mpc_en_m = 3.0856e22   # 1 Mpc =  3.0856e22 m

G = 6.6742*1e-11 * (1/Mpc_en_m)**3 * M_soleil * (Gyr_en_s)**2 # constante gravitationnelle de Newton en Mpc3


Omega_m0 = 0.315
Omega_r0 = 9*1e-5
Omega_L0 = 0.685


mu_x = mu_y = 0
sigma_x = sigma_y = 1




#------------------------------------------------------------------------------------------------------------------------------------
# partie servant a trouver le taux de croissance linéaire D(a)
def f(x):
    return (Omega_m0 * x**-1  + Omega_L0 * x**2  + Omega_r0)**(-3/2)




def H(x):
    return H0 * np.sqrt(Omega_r0 * x**-4  +  Omega_m0 * x**-3  +  Omega_L0)




def integration_rectangle(a, N=1e4):
    Somme = np.zeros( len(a) )

    for i in range( len(a) ):
        if i%10000 == 0:
            print("i =", i)
        x = np.linspace(1e-7, a[i], int(N))
        
        S = 0
        S = np.sum( 0.5 * (x[1:]-x[:-1]) * (f(x[1:]) + 3*f(x[:-1])) )
        
        #for j in range( len(x)-1 ):
        #    S += 0.5 * (x[j+1] - x[j]) * (f(x[j+1]) + f(x[j]))

        Somme[i] = S  
    return Somme   
 




def D(t, h=1e-5):
    x = 1e-5
    I = 0

    A = cpp.a2(t)
    
    while x < A:

        I += 0.5 * h * ( f(x + h) + f(x) ) 
        x += h 
    
    taux_accr = I * H(A) / H0

    
    print(f'D({t}) = ' , taux_accr / 1.000299270341051)   

    return taux_accr / 1.000299270341051         # cste de normalisation : D(a=1) = 1.000299270341051

 

"""
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
    dx = np.abs(x[0,0]-x[0,1])
    dy = np.abs(y[0,0]-y[1,0])



    dpsi_dx = ( psi(x+dx,y, sigma_x, sigma_y, mu_x, mu_y, rho) - psi(x-dx,y, sigma_x, sigma_y, mu_x, mu_y, rho) ) / (2*dx)
    dpsi_dy = ( psi(x,y+dy, sigma_x, sigma_y, mu_x, mu_y, rho) - psi(x,y-dy, sigma_x, sigma_y, mu_x, mu_y, rho) ) / (2*dy)



    
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
"""


    







#------------------------------------------------------------------------------------------------------------------------------------
# Partie simplifiée avec les transformées de Fourier 
def gradient_psi(L_box, N, sc = 0.1):
    dx = dy = L_box/N
    
    phi = np.random.normal(size=(N,N), scale=sc)   # potentiel initial (tableau 'carré'(dim 2) de N points (les valeurs sont aléatoires et reparties selon une loi gaussienne))
    
    psi = phi * -2/(3*Omega_m0*H0**2)

    psi_tilde = np.fft.fftn(psi)  

    kx = 2*np.pi*np.fft.fftfreq( N , dx )
    ky = 2*np.pi*np.fft.fftfreq( N , dy )
    kx,ky = np.meshgrid(kx, ky)
    

    #k = np.sqrt(kx**2 + ky**2)
    #k[0,0] = np.inf
    

    grad_psi_tilde_x = 1j * kx * psi_tilde #/ k
    grad_psi_tilde_y = 1j * ky * psi_tilde #/ k

    grad_psi_x = np.real(np.fft.ifftn(grad_psi_tilde_x))
    grad_psi_y = np.real(np.fft.ifftn(grad_psi_tilde_y))

    return grad_psi_x, grad_psi_y





def position_xy(t, L_box, N):

    qx = np.linspace(-L_box/2, L_box/2, N)
    qy = np.linspace(-L_box/2, L_box/2, N)
    qx,qy = np.meshgrid(qx,qy)
    
    gx, gy = gradient_psi(L_box,N)
    
    px = qx + D(t) * gx
    py = qy + D(t) * gy
    
    return px, py





def affichage_position(t, L_box, N):
    x , y = position_xy(t, L_box, N)
    
    H, xedges, yedges = np.histogram2d(x.ravel(), y.ravel(), bins=200, range=[[-L_box/2, L_box/2],[-L_box/2, L_box/2]])
    
    plt.figure()
    plt.imshow(H, extent=[-L_box/2, L_box/2, -L_box/2, L_box/2], origin='lower', cmap='plasma')
    plt.colorbar(label="Densité")

    if t == ti :
        plt.title(f"Position à t = {t : .6e} Gyr")
    else:
        plt.title(f"Position à t = {t} Gyr")

    plt.xlabel("X en Mpc/h")
    plt.ylabel("Y en Mpc/h")
    
    
    



def affichage_gradient(L_box, N):
    grad_psi_x, grad_psi_y = gradient_psi(L_box, N)
    
    H , xedges, yedges = np.histogram2d(grad_psi_x.ravel(), grad_psi_y.ravel(), bins=200)

    plt.figure()
    plt.imshow(H, extent=[-L_box/2, L_box/2, -L_box/2, L_box/2], origin='lower', cmap='plasma')
    plt.colorbar(label="Densité")
    plt.title("Affichage de grad_psi")
    plt.xlabel("X en Mpc/h")
    plt.ylabel("Y en Mpc/h")
    
    


def plot(T, L_box, N):
    for i in T:
        print('i =', i)
        affichage_position(i, L_box, N)
    plt.show()


L_box = 2000   ### Mpc/h
N = 1500     ### nb de "pixels"


affichage_gradient(L_box, N)
affichage_position(ti, L_box, N)
affichage_position(5, L_box, N)
affichage_position(13.8, L_box, N)
plt.show()





