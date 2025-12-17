import numpy as np
import matplotlib.pyplot as plt
import bib_function as bib



t = np.linspace(1e-7, 14, 500)
L_box = 120  # taille de la boite en Mpc/h
N = 600      # nombre de points dans la grille 



##### courbe du facteur d'échelle en fonction du temps #####
def scalefactor_t():
    plt.figure(figsize=(7,5))
    plt.plot(t, bib.a(t), 'b:')
    plt.xlabel('temps t en Gyr')
    plt.ylabel('facteur d\'échelle a(t)')
    plt.title("Courbe d\'évolution du facteur d\'échelle a(t)")
    plt.grid(True)
    plt.show()





##### courbes des paramètres cosmologiques (omega/rho) #####

# plot des Omega_i(t)
def omega_t():
    plt.figure(figsize=(7,5))
    plt.plot(t, bib.Omega_l(t), 'r:', label="Omega_L(t)")
    plt.plot(t, bib.Omega_m(t), 'b:', label="Omega_M(t)")
    plt.plot(t, bib.Omega_r(t), 'g:', label="Omega_R(t)")
    plt.xlabel("temps t (en Gyr)")
    plt.ylabel("Omega_i(t)")
    plt.title(f"Évolution des paramètres de densité en fonction du temps")
    plt.legend()
    plt.grid(True)
    plt.show()



# plot des Omega_i(a) en log/log: 
def omega_a():
    plt.figure(figsize=(7,5))
    plt.plot(bib.a(t), bib.Omega_l(bib.a(t)), 'r:', label="Omega_L(t)")
    plt.plot(bib.a(t), bib.Omega_m(bib.a(t)), 'b:', label="Omega_M(t)")
    plt.plot(bib.a(t), bib.Omega_r(bib.a(t)), 'g:', label="Omega_R(t)")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("facteur d'échelle a(t)")
    plt.ylabel("Omega_i(t)")
    plt.title(f"Évolution des paramètres de densité en fonction du facteur d'échelle ")
    plt.legend()
    plt.grid(True)
    plt.show()



# plot des rho_i(a) en log/log:
def rho_a():
    plt.figure(figsize=(7,5))
    plt.plot(bib.a(t), bib.Omega_l(bib.a(t))*bib.rho_c(t), 'r:', label="rho_L(t)")
    plt.plot(bib.a(t), bib.Omega_m(bib.a(t))*bib.rho_c(t), 'b:', label="rho_M(t)")
    plt.plot(bib.a(t), bib.Omega_r(bib.a(t))*bib.rho_c(t), 'g:', label="rho_R(t)")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("facteur d'échelle a(t)")
    plt.ylabel("rho_i(t)")
    plt.title(f"Évolution des densités en fonction du facteur d'échelle en log/log")
    plt.legend()
    plt.grid(True)
    plt.show()




# plot des rho_i(t) en log/log:
def rho_t():
    plt.figure(figsize=(7,5))
    plt.plot(t, bib.Omega_l(t)*bib.rho_c(t), 'r:', label="rho_L(t)")
    plt.plot(t, bib.Omega_m(t)*bib.rho_c(t), 'b:', label="rho_M(t)")
    plt.plot(t, bib.Omega_r(t)*bib.rho_c(t), 'g:', label="rho_R(t)")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("temps t en Gyr")
    plt.ylabel("rho_i(t)")
    plt.title(f"Évolution des densités en fonction du temps en log/log")
    plt.legend()
    plt.grid(True)
    plt.show()




##### courbe du facteur d'accroissement #####

# en fonction du temps
def D_t():
    plt.plot(t, bib.D(t))
    plt.xlabel('temps t en Gyr')
    plt.ylabel('growth factor D(t)')
    plt.title('Evolution du growth factor en fonction du temps')
    plt.grid(True)
    plt.show()


# en fonction du facteur d'échelle 
def D_a():
    plt.plot(bib.a(t), bib.D(bib.a(t)))
    plt.xlabel('scale factor a')
    plt.ylabel('growth factor D(a)')
    plt.title('Évolution du growth factor D en fonctin du scale factor a')
    plt.grid(True)
    plt.show()





##### affichage des images #####

# affichage du gradient de psi en 2D
def G_psi():
    bib.affichage_gradient(L_box, N, 0.075)
    plt.show()


#affichage de la paosition en fonction du temps 
def pos_t():
    bib.affichage_position(13.8, L_box, N, 0.5) 
    plt.show()



pos_t()