#include <igl/unproject_onto_mesh.h>
#include <igl/adjacency_list.h>
#include <igl/barycenter.h>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/grad.h>
#include <igl/jet.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/readDMAT.h>
#include <igl/readOFF.h>
#include <igl/repdiag.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <iostream>
#include <sstream>
#include <set>
#include <vector>

#include "../include/imGUIviewer.h"

// using namespace std;
void boundingBox(Eigen::MatrixXd &V_box){}

std::set<int>& get_k_ring_neighbors(std::set<int>& neighbors, const std::vector<std::vector<int>>& adjacency_list,
    int vertex_id, int k)
{
    neighbors.clear();
    std::set<int> current_ring = {vertex_id};
    
    for (int ring = 0; ring < k; ++ring) 
    {
        std::set<int> next_ring;

        for (int v : current_ring) {
            for (int neighbor : adjacency_list[v]) {

                if (neighbors.find(neighbor) == neighbors.end()) {
                    next_ring.insert(neighbor);
                }
            }
        }
        neighbors.insert(current_ring.begin(), current_ring.end());
        current_ring = next_ring;
    }
    neighbors.insert(current_ring.begin(), current_ring.end());
    return neighbors;
}

void anneaux_topologiques(std::set<int>& neighbors, Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXd& neighborhood_colors, int selected_vertex, int k_rings){
        
    std::vector<std::vector<int>> adjacency_list;
    igl::adjacency_list(F, adjacency_list);

    neighbors = get_k_ring_neighbors(neighbors, adjacency_list, selected_vertex, k_rings);
    
    // Visualiser le voisinage
    for (int vid : neighbors) {
        if (vid < V.rows()) {
            neighborhood_colors(vid) = 1.0;
        }
    }
}

bool callback_visualize(igl::opengl::glfw::Viewer& viewer, Eigen::MatrixXd& C, 
    Eigen::MatrixXd& V, Eigen::MatrixXi& F, const Eigen::VectorXd& values){
    
    // Normalise colots
    double min_val = values.minCoeff();
    double max_val = values.maxCoeff();

    double range = max_val - min_val;
    if (std::abs(range) < 1e-9) range = 1.0;
    
    Eigen::MatrixXd colors(V.rows(), 3);
    for (int i = 0; i < V.rows(); ++i) {
        double normalized = (values(i) - min_val) / range;
        // Colormap: bleu -> vert -> rouge
        if (normalized < 0.5) {
            colors(i, 0) = 0.0;                    // Rouge
            colors(i, 1) = 2.0 * normalized;       // Vert
            colors(i, 2) = 1.0 - 2.0 * normalized; // Bleu
        } else {
            colors(i, 0) = 2.0 * (normalized - 0.5); // Rouge
            colors(i, 1) = 1.0 - 2.0 * (normalized - 0.5); // Vert
            colors(i, 2) = 0.0;                     // Bleu
        }
    }
    viewer.data().set_colors(colors);
    std::cout << " (min=" << min_val << ", max=" << max_val << ")" << std::endl;
    return true;
}


bool callback_colorize_face(igl::opengl::glfw::Viewer& viewer, Eigen::MatrixXd& C, 
    Eigen::MatrixXd& V, Eigen::MatrixXi& F){
        int fid;
        Eigen::Vector3f bc;
        // Cast a ray in the view direction starting from the mouse position
        double x = viewer.current_mouse_x;
        double y = viewer.core().viewport(3) - viewer.current_mouse_y;
        if(igl::unproject_onto_mesh(Eigen::Vector2f(x,y), viewer.core().view,
        viewer.core().proj, viewer.core().viewport, V, F, fid, bc))
        {
            // paint hit red
            C.row(fid)<<1,0,0;
            viewer.data().set_colors(C);
            return true;
        }
        return false;
    };

bool callback_select_mode(igl::opengl::glfw::Viewer& viewer, Eigen::MatrixXd& V,
    Eigen::MatrixXi& F, ViewerControls& ui_state)
{
    int fid;
    Eigen::Vector3f bc;

    double x = viewer.current_mouse_x;
    double y = viewer.core().viewport(3) - viewer.current_mouse_y;

    if (igl::unproject_onto_mesh(Eigen::Vector2f(x,y), viewer.core().view,
    viewer.core().proj, viewer.core().viewport, V, F, fid, bc))
    {

        int closest_v_id = -1; //id du vertex le plus proche du clic

        if (bc(0) > bc(1) && bc(0) > bc(2))
        {
            closest_v_id = F(fid, 0); // face d'index fid, sommet 0
        }
        else if (bc(1) > bc(0) && bc(1) > bc(2))
        {
            closest_v_id = F(fid, 1);
        }
        else
        {
            closest_v_id = F(fid, 2);
        }

        ui_state.selected_vertex = closest_v_id;
        viewer.data().clear_points();

        Eigen::MatrixXd P_selected(1, 3);
        P_selected.row(0) = V.row(closest_v_id); //on copie les coord du point dans un vect 1x3

        viewer.data().add_points(P_selected, Eigen::RowVector3d(1, 0, 0)); // couleur rouge

        return true;
    }
    return false;
}

void compute_mean_curvature(Eigen::VectorXd& s_laplacien, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    std::vector<std::vector<int>> adjacency_list;
    igl::adjacency_list(F, adjacency_list);
    s_laplacien.resize(V.rows());

    for (int v = 0; v < V.rows(); ++v)
    {
        int nb_voisin = adjacency_list[v].size();

        if (nb_voisin<1) continue;

        Eigen::RowVector3d moyenne_voisins(0,0,0);
        for (int voisin : adjacency_list[v])
        {
            moyenne_voisins += V.row(voisin);
        }
        moyenne_voisins /= nb_voisin;
        
        Eigen::RowVector3d vec_Laplacien = V.row(v) - moyenne_voisins;
        s_laplacien(v) = vec_Laplacien.norm();
    }
}



void step_heat_diffusion(Eigen::VectorXd& U, const std::vector<std::vector<int>>& adj_list, std::set<int>& neighbors)
{
    Eigen::VectorXd U_next = U;

    bool use_neighbors = !neighbors.empty();

    for (int v = 0; v < U.rows(); ++v) 
    {
        if (use_neighbors && neighbors.find(v) == neighbors.end()) continue;

        std::vector<int> voisins = adj_list[v];

        if (voisins.empty()) continue;
    
        double moyenne = 0; 
        for (int voisin : adj_list[v])
        {
            moyenne += U(voisin); 
        }
        moyenne /= adj_list[v].size();

        double laplacien = moyenne - U(v); 

        U_next(v) = U(v) + 0.5 * laplacien; 
    }
    U = U_next;
}



int main(int argc, char *argv[])
{
    // V for vertices and F for faces
    Eigen::MatrixXd V, C, U;
    Eigen::MatrixXi F;
    Eigen::SparseMatrix<double> L;
    
    /// Load Model //
    // viewer.load_mesh_from_file(".../data/sphere.obj");

    // Load a mesh in OFF format //
    if (!igl::readOFF("../data/cow.off", V, F)) {
        std::cerr << "Erreur: impossible de charger le maillage" << std::endl;
        return -1;
    }

    igl::cotmatrix(V,F,L);

    //////////////////////////
    ///// Init program ///////
    //////////////////////////

    std::vector<std::vector<int>> adj_list;
    igl::adjacency_list(F, adj_list);

    Eigen::VectorXd heat_values = Eigen::VectorXd::Zero(V.rows());
    
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
    
    std::set<int> neighbors;

    Eigen::VectorXd s_laplacien(V.rows());

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

    Eigen::VectorXd neighborhood_colors = Eigen::VectorXd::Zero(V.rows());
    
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
                }

                for(int i=0; i< 10; i++){ // pour aller plus vit
                
                    step_heat_diffusion(heat_values, adj_list, neighbors);
                    callback_visualize(viewer, C, V, F, heat_values);
                }
                std::cout << "Diffusion en cours..." << std::endl;

                // else {
                //     for(int i=0; i< 10; i++){
                        
                //         step_heat_diffusion(heat_values, adj_list);       
                //         callback_visualize(viewer, C, V, F, heat_values);
                //     }   
                // }                
                return true;
            }
            case '1':
                ui_state.selection_mode = true; 
                viewer.data().set_colors(C);
                viewer.data().clear_points(); 
                std::cout << "Mode Selection ACTIVE. Cliquez sur le maillage." << std::endl;
                return true;

            // case '1': //selection face
            //     viewer.data().set_colors(C);
            //     viewer.data().clear_points(); // on efface le pt quand on change de mode
            //     return true;

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
                compute_mean_curvature(s_laplacien, V, F) ;

                // Eigen::SparseMatrix<double> M; // Matrice de masse
                
                // igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
                
                // Eigen::SparseMatrix<double> Minv;
                // Minv = M.cwiseInverse(); 
                // Eigen::SparseMatrix<double> L_normalized = Minv * L;
                
                // // Calculer les valeurs du Laplacien pour chaque sommet
                // Eigen::VectorXd s_laplacien = L_normalized * V.col(0); 
                
                callback_visualize(viewer, C, V, F, s_laplacien);
                return true;
            }

            case '4': //lissage laplacien
            {
                Eigen::SparseMatrix<double> M;
                igl::massmatrix(U,F,igl::MASSMATRIX_TYPE_BARYCENTRIC,M);
                
                // Solve (M-delta*L) U = M*U
                const auto & S = (M - 0.001*L);
                
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
                
                // break;
                
                // Send new positions, update normals, recenter
                viewer.data().set_vertices(U);
                viewer.data().compute_normals();
                viewer.core().align_camera_center(U,F);
                return true;
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
    // viewer.data().add_label(viewer.data().V.row(0) + viewer.data().V_normals.row(0).normalized()*0.005, "Hello World!");
    viewer.launch();
    
    return 0;
}
