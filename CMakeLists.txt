# CMakeLists.txt pour un projet utilisant libigl
# Version minimale moderne et complète

cmake_minimum_required(VERSION 3.16)

# Nom du projet
project(mesh_tp)

# Configuration C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Configuration des types de build
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Configuration des flags de compilation
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# =====================================
# MÉTHODE 1: Téléchargement automatique (RECOMMANDÉ)
# =====================================

include(FetchContent)

# Télécharger libigl automatiquement
FetchContent_Declare(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG v2.5.0
)

# Rendre libigl disponible
FetchContent_MakeAvailable(libigl)

# =====================================
# MÉTHODE 2: Si vous avez libigl localement (ALTERNATIVE)
# =====================================
# Décommentez les lignes suivantes si vous préférez une copie locale:
# set(LIBIGL_ROOT "/chemin/vers/votre/libigl")
# add_subdirectory(${LIBIGL_ROOT} libigl)

# =====================================
# Configuration des modules libigl nécessaires
# =====================================

# Module pour visualisation (Viewer avec GLFW)
igl_include(core)      # Fonctionnalités de base
igl_include(opengl)    # Support OpenGL (requis pour ImGui)
igl_include(glfw)
igl_include(imgui)

# Module pour les opérations géométriques avancées (optionnel)
# igl_include(copyleft cgal)

# Module pour les maillages triangulaires (optionnel)
# igl_include(restricted triangle)

# =====================================
# Création de l'exécutable
# =====================================

# Fichier source principal
add_executable(${PROJECT_NAME} main.cpp)

# =====================================
# Linkage des bibliothèques
# =====================================

# Lier avec les modules libigl nécessaires
target_link_libraries(${PROJECT_NAME} PUBLIC 
    igl::core
    igl::opengl  
    igl::glfw    
    igl::imgui 
)

# Si vous utilisez des modules optionnels, ajoutez-les:
# target_link_libraries(${PROJECT_NAME} PUBLIC igl_copyleft::cgal)
# target_link_libraries(${PROJECT_NAME} PUBLIC igl_restricted::triangle)

# =====================================
# Configuration additionnelle (optionnel)
# =====================================

# Dossier pour les ressources (maillages, textures, etc.)
file(COPY ${CMAKE_SOURCE_DIR}/data DESTINATION ${CMAKE_BINARY_DIR})

# Définir le chemin vers les données
target_compile_definitions(${PROJECT_NAME} PRIVATE 
    DATA_PATH="${CMAKE_BINARY_DIR}/data/"
)

# =====================================
# Configuration de débogage (optionnel)
# =====================================

# Afficher les informations de compilation
message(STATUS "Projet: ${PROJECT_NAME}")
message(STATUS "Type de build: ${CMAKE_BUILD_TYPE}")
message(STATUS "Compilateur C++: ${CMAKE_CXX_COMPILER}")
message(STATUS "Version C++: ${CMAKE_CXX_STANDARD}")

# =====================================
# Pour votre TP spécifiquement
# =====================================

# Exemple de fichier source pour votre TP
# Remplacez main.cpp par votre fichier principal

# Si vous avez plusieurs fichiers sources:
# set(SOURCES 
#     main.cpp
#     mesh_processing.cpp
#     laplacian_computation.cpp
# )
# add_executable(${PROJECT_NAME} ${SOURCES})