"""

Bibliothèque de fonctions pour le projet :
- facteur d’échelle a(t)
- paramètres de densité
- taux de croissance D(t)
- approximation de Zel’dovich
- affichage

"""




import numpy as np 
import matplotlib.pyplot as plt



# variables globales
H0 = 0.070                          # 67.4*1e-19*1e9*365.25*24*3600/3.03 : cste de Hubble à notre ère (conversion en Gyr-1)
ai = 1e-6                           # facteur d'échelle initial (au début de l'Univers) 
ti = 1/(365.25*1e9)                 # i = temps initial (vers les débuts de l'Univers) en Gyr
t0 =  14                            # 0 = temps de notre ère (age de l'Univers) en Gyr
                   
Omega_m0 = 0.315
Omega_r0 = 9*1e-5
Omega_L0 = 0.685

M_soleil = 1.98847e30  # masse du soleil en kg
Gyr_en_s = 3.15576e16  # 1 Gyr = 3.15576e16 s
Mpc_en_m = 3.0856e22   # 1 Mpc =  3.0856e22 m

G = 6.6742 * 1e-11 * (1/Mpc_en_m)**3 * M_soleil * (Gyr_en_s)**2     # constante gravitationnelle de Newton en Mpc3





# facteur d'échelle a(t) :

def a(t_voulu, om_0=Omega_m0, or_0=Omega_r0, ol_0=Omega_L0, H = 1e-5, afficher = False):
    norm = 1.0154729584922937    # cste de normalisation : a_non_normé(t0) = 1.0154729584922937
    A = ai
    t = ti

    if np.isscalar(t_voulu):

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
            print(f'a({t_voulu}) =', A/norm)

        return A/norm
    
    else : 
        L = np.zeros( len(t_voulu) )
        A = ai
        t = ti

        for j, tj in enumerate(t_voulu):
            

            while t < tj:

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

            L[j] = A
        
        return L/norm







######## paramètres cosmologiques ########

# paramètre de densité de la matiere totale (mat visible + mat noire)
def Omega_m(t, om_0=Omega_m0, or_0=Omega_r0, ol_0=Omega_L0):  

    Omega_m = om_0 / ( (a(t)**-1)*or_0 + om_0 + (a(t)**3) * ol_0 )
    
    return Omega_m



# paramètre de densité de radiation
def Omega_r(t, om_0=Omega_m0, or_0=Omega_r0, ol_0=Omega_L0): 

    Omega_r = or_0 / ( or_0 + a(t) * om_0 + (a(t)**4) * ol_0 ) 
    
    return Omega_r



# paramètre de densité d'énergie noire 
def Omega_l(t, om_0=Omega_m0, or_0=Omega_r0, ol_0=Omega_L0):     

    Omega_l = ol_0 / ( (a(t)**-4) * or_0 + (a(t)**-3) * om_0 + ol_0 ) 
    
    return Omega_l



# densité critique
def rho_c(t, om_0=Omega_m0, or_0=Omega_r0, ol_0=Omega_L0):

    rho_c = (3/(8*np.pi*G)) * H0**2 * ((a(t)**-4) * or_0 + (a(t)**-3) * om_0 + ol_0)
    
    return  rho_c




####### détermination du taux de croissance linéaire des structures D(t) ########

# intégrande dans l'eq de D(t)
def f(x):
    """Intégrande pour le calcul de D(a)"""
    return (Omega_m0 * x**-1  + Omega_L0 * x**2  + Omega_r0)**(-3/2)



# équation de Friedmann
def H(x):
    """Équation de Friedmann : H(a)"""
    return H0 * np.sqrt(Omega_r0 * x**-4  +  Omega_m0 * x**-3  +  Omega_L0)



# taux de croissance linéaire des structures D(t)
def D(t, h=1e-5):

    """
    Calcule le taux de croissance linéaire D(t) via l'intégrale :   D(a) ∝ H(a) ∫ dx / (x^3 H(x)^3)
    
    Normalisation imposée pour D(a=1) = 1.

    Gère scalaires et tableaux.
    """

    norm = 1.000299270341051                           # cste de normalisation : D_non_normalisée(a=1) = 1.000299270341051
    
    # cas scalaire
    if np.isscalar(t) :   
        x = 1e-5
        I = 0
        A = a(t)
        
        while x < A:

            I += 0.5 * h * ( f(x + h) + f(x) )        # méthode des trapèzes
            x += h 
    
        taux_accr = I * H(A) / H0

        print(f'D({t}) = ' , taux_accr / norm)   

        return taux_accr / norm   
          

    # cas tableau
    else :             
        L = np.zeros( len(t) )
        
        for i, ti in enumerate(t): 
            Ai = a(ti)
            x = 1e-5
            I = 0
            
            while x < Ai:

                I += 0.5 * h * ( f(x + h) + f(x) ) 
                x += h 
    
            taux_accr = I * H(Ai) / H0   

            L[i] = taux_accr / norm         

        return L



# calcul du gradient de psi :
def gradient_psi(L_box, N, sc = 1):

    """
    Calcul du gradient de psi en utilisant des transformées de Fourier
    - phi : potentiel initial aléatoire (gaussien)
    - psi : potentiel gravitationnel lié par la relation de Poisson en espace de Fourier
    - la dérivée est obtenue via multiplication par i*k
    """

    dx = dy = L_box/N                                   # maillage carré de taille L_box et de N x N points

    kx = 2*np.pi*np.fft.fftfreq( N , dx )               # vecteurs d'ondes kx et ky correspondant aux fréquences de la FFT ( k_x = 2*pi/lambda )
    ky = 2*np.pi*np.fft.fftfreq( N , dy )
    kx,ky = np.meshgrid(kx, ky)
   

    k = np.sqrt(kx**2 + ky**2)
    k[0,0] = np.inf                                     # pour éviter une division par zéro quand k=0

    phi = np.random.normal(size=(N,N), scale=sc)        # potentiel initial (tableau 'carré'(dim 2) de NxN points (les valeurs sont aléatoires et reparties selon une loi gaussienne avec écart-type "sc" )
    phi_tilde = np.fft.fftn(phi) 


    psi_tilde = phi_tilde * -2/(3* k * Omega_m0*H0**2)


    grad_psi_tilde_x = 1j * kx * psi_tilde              # dérivée dans l'espace de Fourier
    grad_psi_tilde_y = 1j * ky * psi_tilde

    grad_psi_x = np.real(np.fft.ifftn(grad_psi_tilde_x))
    grad_psi_y = np.real(np.fft.ifftn(grad_psi_tilde_y))

    return grad_psi_x, grad_psi_y





# détermination de la position au cours du temps :
def position_xy(t, L_box, N, sc):

    """
    Approximation de Zel’dovich :   x(t) = q + D(t) * nabla_psi(q)
    
    q : grille régulière (coordonnées initiales)
    """

    qx = np.linspace(-L_box/2, L_box/2, N)
    qy = np.linspace(-L_box/2, L_box/2, N)
    qx,qy = np.meshgrid(qx,qy)
    
    gx, gy = gradient_psi(L_box,N, sc)
    
    px = qx + D(t) * gx
    py = qy + D(t) * gy
    
    return px, py





######## affichage des graphes ########

# représentation graphique de la position au cours du temps t :
def affichage_position(t, L_box, N, sc):
    x , y = position_xy(t, L_box, N, sc)
    
    H, xedges, yedges = np.histogram2d(x.ravel(), y.ravel(), bins=300, range=[[-L_box/2, L_box/2],[-L_box/2, L_box/2]])     # carte de densité
    
    plt.figure()
    plt.imshow(H, extent=[-L_box/2, L_box/2, -L_box/2, L_box/2], origin='lower', cmap='inferno')
    plt.colorbar(label="Densité")

    if t == ti :
        plt.title(f"Position à t = {t : .6e} Gyr, N_px = {N},  scale = {sc}")
    else:
        plt.title(f"Position à t = {t} Gyr, N_px = {N}, scale = {sc}")

    plt.xlabel("X en Mpc/h")
    plt.ylabel("Y en Mpc/h")
    
    

# représentation graphique du gradient de psi
def affichage_gradient(L_box, N, sc):
    grad_psi_x, grad_psi_y = gradient_psi(L_box, N, sc)
    
    H , xedges, yedges = np.histogram2d(grad_psi_x.ravel(), grad_psi_y.ravel(), bins=250)

    plt.figure()
    plt.imshow(H, extent=[-L_box/2, L_box/2, -L_box/2, L_box/2], origin='lower', cmap='inferno')
    plt.colorbar(label="Densité")
    plt.title(f"Affichage de grad_psi , scale = {sc}")
    plt.xlabel("X en Mpc/h")
    plt.ylabel("Y en Mpc/h")
    