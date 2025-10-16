#include <igl/readOFF.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/adjacency_list.h>
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
std::set<int> get_k_ring_neighbors(const std::vector<std::vector<int>>& adjacency_list, int vertex_id, int k) {
    std::set<int> neighbors;
    std::set<int> current_ring = {vertex_id};
    
    for (int ring = 0; ring < k; ++ring) {
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
bool callback_click(igl::opengl::glfw::Viewer& viewer, Eigen::MatrixXd& C, 
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

int main(int argc, char *argv[])
{
    // V for vertices and F for faces
    Eigen::MatrixXd V, C;
    Eigen::MatrixXi F;
    
    /// Load Model
    // viewer.load_mesh_from_file(".../data/sphere.obj");
    // Load a mesh in OFF format
    if (!igl::readOFF("../data/bunny.off", V, F)) {
        std::cerr << "Erreur: impossible de charger le maillage" << std::endl;
        return -1;
    }
    
    ////////////////////////////
    ///// COLOR SETTINGS //////
    ///////////////////////////
    C = Eigen::MatrixXd::Constant(F.rows(),3,1);

    Eigen::Vector3d m = V.colwise().minCoeff();
    Eigen::Vector3d M = V.colwise().maxCoeff();

    
    // Affichage //
    igl::opengl::glfw::Viewer viewer;
    
    
    //////////////////
    /// CALLBACK /////
    //////////////////
    
    // Define the callback for face capture
    viewer.callback_mouse_down = [&](igl::opengl::glfw::Viewer& viewer, int button, int modifier) {
        return callback_click(viewer,C , V, F);
    };
    
    
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

    /////////////////////////
    ////// TP et Projet /////
    /////////////////////////
    
    double doubleVariable = 0.1f;


    /////////////////////////////////////////////////////
    // Selection du voisinage par anneaux topologiques //
    /////////////////////////////////////////////////////

    std::vector<std::vector<int>> adjacency_list;
    igl::adjacency_list(F, adjacency_list);
    int selected_vertex = 0;
    int k_rings = 2;
    auto neighbors = get_k_ring_neighbors(adjacency_list, selected_vertex, k_rings);

    /////////////////////////
    ///////// PLOT //////////
    /////////////////////////

    viewer.data().set_mesh(V, F);
    viewer.data().set_colors(C);
    // viewer.data().add_label(viewer.data().V.row(0) + viewer.data().V_normals.row(0).normalized()*0.005, "Hello World!");
    viewer.launch();
}
