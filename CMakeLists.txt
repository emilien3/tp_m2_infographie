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
# MÉTHODE 1: Téléchargement automatique
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
# Configuration des modules libigl nécessaires
# =====================================
# Module pour visualisation (Viewer avec GLFW)
igl_include(core)      # Fonctionnalités de base
igl_include(opengl)    # Support OpenGL (requis pour ImGui)
igl_include(glfw)
igl_include(imgui)

include_directories(${CMAKE_SOURCE_DIR}/include)

# =====================================
# Création de l'exécutable
# =====================================

file( GLOB SOURCES
    "src/*.cpp"
)
add_executable(${PROJECT_NAME} ${SOURCES})



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

# =====================================
# Configuration additionnelle (optionnel)
# =====================================
# Dossier pour les ressources (maillages, textures, etc.)
file(COPY ${CMAKE_SOURCE_DIR}/data DESTINATION ${CMAKE_BINARY_DIR})

# Définir le chemin vers les données
target_compile_definitions(${PROJECT_NAME} PRIVATE 
    DATA_PATH="${CMAKE_BINARY_DIR}/data/"
)


message(STATUS "Projet: ${PROJECT_NAME}")
message(STATUS "Type de build: ${CMAKE_BUILD_TYPE}")
message(STATUS "Compilateur C++: ${CMAKE_CXX_COMPILER}")
message(STATUS "Version C++: ${CMAKE_CXX_STANDARD}")

