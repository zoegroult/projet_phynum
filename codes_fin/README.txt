
1. Présentation générale 
-------------------------

Ce projet a pour objectif de modéliser l’évolution de l’Univers depuis la recombinaison (CMB) jusqu'à la formation des grandes structures (cosmic web). 
Le code permet de calculer le facteur d’échelle, les paramètres de densité cosmologiques, la densité critique ainsi que le taux de croissance linéaire des structures.
Il inclut également une simulation simplifiée de la formation des structures en 2D via l’approximation de Zel’dovich et des outils de visualisation associés.




2. Description des fichiers 
----------------------------

- bib_function.py :

  Bibliothèque contenant l’ensemble des fonctions necessaires au projet :
    • facteur d’échelle a(t) ,
    • paramètres de densité Ω_m, Ω_r, Ω_Λ ,
    • densité critique ρ_c , 
    • facteur de croissance D(t) ,
    • calcul du potentiel gravitationnel et approximation de Zel’dovich.

  Ce fichier ne produit pas de résultats graphiques par lui-même.


- code_graphes.py :

  Script de visualisation du projet.
  Il regroupe les fonctions permettant de tracer les différentes grandeurs physiques
  (facteur d’échelle, paramètres de densité, densités, facteur de croissance)
  ainsi que les cartes de densité issues de l’approximation de Zel’dovich.
  Il utilise les fonctions définies dans bib_function.py.





3. Instructions d’exécution
----------------------------

Prérequis :
    - Python 3
    - Bibliothèques : numpy, matplotlib
    - Fichiers : bib_function.py et code_graphes.py


Exécution :
    - Lancer code_graphes.py
    - Un menu interactif s’affichera, permettant de choisir quel plot visualiser. 
        Le menu interactif (fonctionnalité ajoutée avec l’aide de ChatGPT) a été conçu pour faciliter la lecture et l’utilisation du 
        code pour les correcteurs, afin qu’ils puissent visualiser un plot à la fois sans avoir tous les graphes affichés en même temps.  
        Voici son fonctionnement :
            • Entrer le numéro correspondant au plot souhaité dans le menu.
            • Fermer la figure pour revenir au menu.
            • Entrer 0 pour quitter le programme.

Personnalisation
Les paramètres (temps, taille de la boîte, nombre de points, scale, etc.) sont définis au début du script et peuvent être modifiés si nécessaire.


