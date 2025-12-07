
# =====================================
# UTILISATION
# =====================================

0. git clone https://github.com/libigl/libigl.git
1. Créez un dossier build: mkdir build && cd build
2. Configurez: cmake ..
3. Compilez: make -j4 (Linux/Mac) ou cmake --build . (multiplateforme)
4. Exécutez: ./mesh_tp

# =====================================


Utilisation du code :


* Cliquer sur 1 : Active le mode de sélection
* Cliquer sur 0 : permet de reset les modifications apporté aux couleurs ou au maillage
* Cliquer avec la souris permet de sélectionner un point

Après ces étapes plusieurs choix (cliquer sur 0 puis avec la souris pour tester d'autres fonctionnalités) :

* Sélectionner dans le bandeau les k rings et cliquer sur 2 pour les afficher
* Cliquer sur 3 pour afficher la courbure moyenne sur l'objet 
* Cliquer sur 4 à plusieurs reprise pour affiner l'objet 

* Pour l'affichage du laplacien :

  * Topologique : 
    * Explicite : Barre Espace (avec ou sans krings)
    * Implicite : Appuie sur 5 pour un saut

  * Géométrique :
    * Diffusion : Appuie sur 6 pour une étape
    * Par résolution de systèmes d'équation : Appuie sur 7
  
* Pour la déformation : 
  * Appuie sur 8
  * Sélection dans la bar de la fonction de transfert  