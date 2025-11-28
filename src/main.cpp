#include <igl/adjacency_list.h>
#include <igl/barycenter.h>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>

// Pas utilisé dans le code 
#include <igl/per_vertex_normals.h>
#include <igl/grad.h>
#include <igl/jet.h>
#include <igl/readDMAT.h>
#include <igl/repdiag.h>
#include <igl/slice_into.h>


#include <igl/massmatrix.h>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include <Eigen/Sparse>

#include "../include/imGUIviewer.h"
#include "../include/utils.hpp"


int main(int argc, char *argv[])
{
    // V for vertices and F for faces
    Eigen::MatrixXd V, C, U;
    Eigen::MatrixXi F;

    Eigen::SparseMatrix<double> P, L_cot, Mass_bary;
    Eigen::SparseMatrix<double> M_voronoi; // Matrice de masse

    std::vector<std::vector<int>> adj_list;
    std::set<int> neighbors;
    
    /// Load Model //
    // viewer.load_mesh_from_file(".../data/sphere.obj");
    // Load a mesh in OFF format //
    if (!igl::readOFF("../data/cow.off", V, F)) {
        std::cerr << "Erreur: impossible de charger le maillage" << std::endl;
        return -1;
    }
    
    //////////////////////////
    ///// Init program ///////
    //////////////////////////

    igl::adjacency_list(F, adj_list);

    Eigen::VectorXd heat_values = Eigen::VectorXd::Zero(V.rows());    
    
    igl::cotmatrix(V, F, L_cot);
    // igl::cotmatrix(V,F,L_cot);
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, Mass_bary);
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M_voronoi);


                
    
    Eigen::VectorXd s_laplacien(V.rows());
    Eigen::VectorXd neighborhood_colors = Eigen::VectorXd::Zero(V.rows());
    
    ////////////////////////////
    ///// COLOR SETTINGS //////
    ///////////////////////////
    C = Eigen::MatrixXd::Constant(F.rows(),3,1);
    Eigen::Vector3d m = V.colwise().minCoeff();
    Eigen::Vector3d M = V.colwise().maxCoeff();

    // Affichage //
    igl::opengl::glfw::Viewer viewer;
    
    /////////////////
    ///// ImGUI /////
    ////////////////
    
    // Rendering of text labels is handled by ImGui, so we need to enable the ImGui
    // plugin to show text labels.
    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    viewer.plugins.push_back(&plugin);

    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);
    
    ViewerControls ui_state;
    menu.callback_draw_viewer_window = [&]() {
        ui_state.render(menu);
    };

    //////////////////
    /// CALLBACKS ////
    //////////////////
    
    // Define the callback for face capture
    viewer.callback_mouse_down = 
    [&](igl::opengl::glfw::Viewer& viewer, int button, int modifier) -> bool
    {    
        if (ui_state.selection_mode)
        {
            
            return callback_select_mode(viewer, V, F, ui_state);

        }
        // else
        // {
        //     return callback_colorize_face(viewer, C, V, F);
        // }
    };

    viewer.callback_key_down =
    [&](igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) -> bool
    {

        std::cout << "Le mode est à" << ui_state.selection_mode;

        switch (key)
        {
            case '0': // on reset
            {
                if (ui_state.selected_vertex == -1) {
                    std::cout << "Erreur: Selectionnez un sommet d'abord (clic droit) !" << std::endl;
                    return true;
                }
                
                heat_values.setZero(); 
                
                heat_values(ui_state.selected_vertex) = 1.0; 
                
                std::cout << "Source de chaleur placee sur le sommet " << ui_state.selected_vertex << std::endl;
                U = V;

                // Send new positions, update normals, recenter
                viewer.data().set_vertices(U);
                viewer.data().compute_normals();
                viewer.core().align_camera_center(U,F);

                // On affiche
                callback_visualize(viewer, C, V, F, heat_values);
                return true;
            }

            case ' ': // Diffusion
            {

                if (ui_state.boundary){
                    int k_rings = ui_state.k_rings;
                    int selected_vertex = ui_state.selected_vertex;
                    neighbors = get_k_ring_neighbors(neighbors, adj_list, selected_vertex, k_rings);
                }else {
                    neighbors.clear(); // Si boundary désactivé, on diffuse partout
                }

                for(int i=0; i< 10; i++){ // pour aller plus vit
                    step_heat_diffusion(heat_values, adj_list, neighbors);

                    if (ui_state.selected_vertex != -1) heat_values(ui_state.selected_vertex) = 1.0;
                }
                callback_visualize(viewer, C, V, F, heat_values);
                std::cout << "Diffusion en cours..." << std::endl;
                                
                return true;
            }
            case '1':
                ui_state.selection_mode = true; 
                viewer.data().set_colors(C);
                viewer.data().clear_points(); 
                std::cout << "Mode Selection ACTIVE. Cliquez sur le maillage." << std::endl;
                return true;

            case '2':  // k-rings
            {              
                Eigen::VectorXd neighborhood_colors = Eigen::VectorXd::Zero(V.rows());
                int k_rings = ui_state.k_rings;
                int selected_vertex = ui_state.selected_vertex;
                anneaux_topologiques(neighbors, V, F, neighborhood_colors, selected_vertex, k_rings);
                
                std::cout << "Affichage: \n Voisinage " + std::to_string(k_rings) + "-rings";
                callback_visualize(viewer, C, V, F, neighborhood_colors);
                return true;
            }
            
            case '3': //laplacien courbure moyenne
            {
                // compute_mean_curvature(s_laplacien, V, F) ;
                
                Eigen::SparseMatrix<double> Minv;
                Minv = M_voronoi.cwiseInverse(); 
                Eigen::SparseMatrix<double> L_normalized = Minv * L_cot;
                
                // Calculer les valeurs du Laplacien pour chaque sommet
                Eigen::MatrixXd v_norm = L_normalized * V;
                
                Eigen::VectorXd s_laplacien = v_norm.rowwise().norm();//norme ligne par ligne pour l'affichage
                
                callback_visualize(viewer, C, V, F, s_laplacien);
                return true;
            }

            case '4': //lissage laplacien
            {
                Eigen::SparseMatrix<double> M;
                igl::massmatrix(U,F,igl::MASSMATRIX_TYPE_BARYCENTRIC,M);
                
                // Solve (M-delta*L_cot) U = M*U
                const auto & S = (M - 0.001*L_cot);
                
                Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver(S);
                assert(solver.info() == Eigen::Success);
                U = solver.solve(M*U).eval();

                // Compute centroid and subtract (also important for numerics)
                Eigen::VectorXd dblA;
                igl::doublearea(U,F,dblA);
                double area = 0.5*dblA.sum();
                Eigen::MatrixXd BC;
                igl::barycenter(U,F,BC);
                Eigen::RowVector3d centroid(0,0,0);
                for(int i = 0;i<BC.rows();i++)
                {
                centroid += 0.5*dblA(i)/area*BC.row(i);
                }
                U.rowwise() -= centroid;

                // Normalize to unit surface area (important for numerics)
                U.array() /= sqrt(area);
                                
                // Send new positions, update normals, recenter
                viewer.data().set_vertices(U);
                viewer.data().compute_normals();
                viewer.core().align_camera_center(U,F);
                return true;
            }

            case '5': // laplacien par sytèmes d'équation 
            {    

                if (ui_state.boundary){
                    int k_rings = ui_state.k_rings;
                    int selected_vertex = ui_state.selected_vertex;
                    neighbors = get_k_ring_neighbors(neighbors, adj_list, selected_vertex, k_rings);
                } else {
                    neighbors.clear(); // Si boundary désactivé, on diffuse partout
                }
                double lambda = 10; 

                if (ui_state.boundary){
                    int k_rings = ui_state.k_rings;
                    int selected_vertex = ui_state.selected_vertex;
                    neighbors = get_k_ring_neighbors(neighbors, adj_list, selected_vertex, k_rings);
                }else {
                    neighbors.clear(); // Si boundary désactivé, on diffuse partout
                }
                
                
                build_laplacian_matrix(V, F, P, neighbors, ui_state.selected_vertex);
                
                Eigen::VectorXd U_new;
                std::cout << "Diffusion implicite (lambda=" << lambda << ")..." << std::endl;
                
                if (ui_state.selected_vertex != -1) heat_values(ui_state.selected_vertex) = 1.0;
                solve_implicit_heat(heat_values, P, lambda, U_new);
                
                heat_values = U_new;
                callback_visualize(viewer, C, V, F, heat_values);
                return true;
            }

            case '6': // DIFFUSION COTANGENT
            {
                double lambda = 10000.0; 

                Eigen::VectorXd U_new = heat_values;
                
                std::cout << "Diffusion Cotangent (Physique)..." << std::endl;
                
                if (ui_state.selected_vertex != -1) heat_values(ui_state.selected_vertex) = 1.0;

                // Résolution
                solve_cotangent_heat(heat_values, L_cot, Mass_bary, lambda, U_new);
                
                heat_values = U_new;
                
                if (ui_state.selected_vertex != -1) heat_values(ui_state.selected_vertex) = 1.0;

                callback_visualize(viewer, C, V, F, heat_values);
                return true;
            }
            
            case '7':
            {
                if (ui_state.selected_vertex == -1) {
                    std::cout << "Erreur: Selectionnez un point source ('1' -> Clic) !" << std::endl;
                    return true;
                }

                if (ui_state.boundary){
                    neighbors = get_k_ring_neighbors(neighbors, adj_list, ui_state.selected_vertex, ui_state.k_rings);
                } else {
                    neighbors.clear(); 
                }

                Eigen::SparseMatrix<double> L_in_in, L_in_b;
                
                Eigen::VectorXi all, in, b, IA, IC;
                
                Eigen::VectorXd bc, U_in;
                
                solve_laplace_system(L_cot, V, F, L_in_in, L_in_b, all, in, b, IA, IC, bc, ui_state, neighbors, U_in);

                igl::slice_into(U_in, in, heat_values);
                igl::slice_into(bc, b, heat_values);
                callback_visualize(viewer, C, V, F, heat_values);
            }

            default :
                std::cout << "Key is not defined yet";
                return false; 
        }
        
    };
    
    /////////////////////////
    ///////// PLOT //////////
    /////////////////////////

    U = V;
    viewer.data().set_mesh(U, F);
    viewer.data().set_colors(C);
    viewer.launch();
    
    return 0;
}
