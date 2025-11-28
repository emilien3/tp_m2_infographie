
---

### Axe 1 : La Qualité de la Matrice (Le "Moteur")
*Comment on calcule les poids entre les voisins ?*

1.  **Le Laplacien Uniforme (Ce que tu as codé)**
    * **Formule :** "J'ai 5 voisins, je prends $1/5$ de chacun".
    * **Avantage :** Très facile à comprendre et à coder soi-même.
    * **Inconvénient :** Pas très précis physiquement (ne voit pas si un triangle est grand ou petit).
    * **Utilisation :** TP, tests rapides, lissage simple.

2.  **Le Laplacien Cotangent (Celui de LibIGL / Case '4')**
    * **Formule :** Utilise les angles et les aires des triangles.
    * **Avantage :** C'est la vérité physique parfaite.
    * **Inconvénient :** Complexe à coder (besoin de librairies).
    * **Utilisation :** Simulations physiques sérieuses, recherche.

---

### Axe 2 : La Méthode de Résolution (La "Pédale d'accélération")
*Comment on avance dans le temps ?*

1.  **Explicite (Touche ESPACE)**
    * **Méthode :** Boucle `for`. On calcule l'état $t+1$ à partir de $t$.
    * **Analogie :** Marcher pas à pas.
    * **Problème :** Si tu fais un trop grand pas, tu tombes (le système explose/crash).
    * **Code :** `U = U + lambda * L * U`

2.  **Implicite (Touche '5')**
    * **Méthode :** Système Linéaire ($Ax=b$). On demande à l'ordi de trouver l'état d'équilibre futur.
    * **Analogie :** Se téléporter à la destination.
    * **Avantage :** Tu ne tombes jamais. Tu peux faire des sauts de géant (Lambda énorme).
    * **Code :** `solver.solve(...)`

---

### Axe 3 : Le But Physique (La "Destination")
*Qu'est-ce qu'on cherche à obtenir ?*

1.  **Diffusion (Chaleur)**
    * **Équation :** $\frac{\partial u}{\partial t} = \Delta u$
    * **Visuel :** Une couleur qui bave et s'étale.
    * **Data :** On modifie un vecteur de couleurs `Heat`. La forme du lapin ne bouge pas.

2.  **Lissage (Géométrie)**
    * **Équation :** $\frac{\partial P}{\partial t} = \Delta P$
    * **Visuel :** Le lapin fond, les oreilles disparaissent.
    * **Data :** On modifie la matrice `V` (positions).

3.  **Harmonique (Équilibre)**
    * **Équation :** $\Delta u = 0$ (avec bords fixés).
    * **Visuel :** Un dégradé parfait et statique (comme une bulle de savon). Pas de mouvement.
    * **Data :** On cherche l'état final direct pour remplir un trou ou interpoler.

---

### Résumé de ton TP pour t'y retrouver



Actuellement, ton code contient :

* **Touche ESPACE :**
    * Moteur : **Uniforme** (Ta fonction).
    * Méthode : **Explicite** (Boucle).
    * But : **Diffusion** de chaleur.

* **Touche '5' :**
    * Moteur : **Uniforme** (Ta fonction, via la matrice `P`).
    * Méthode : **Implicite** (Solver Eigen).
    * But : **Diffusion** de chaleur.

* **Touche '4' (Code copié-collé) :**
    * Moteur : **Cotangent** (`igl::cotmatrix`).
    * Méthode : **Implicite** (`SimplicialLLT`).
    * But : **Lissage** géométrique (modif `V`).

Tu as simplement exploré toutes les combinaisons possibles de ce "couteau suisse" mathématique ! C'est normal d'être un peu perdu, c'est la richesse du sujet.