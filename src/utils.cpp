#include <igl/adjacency_list.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>

#include <igl/slice.h>
#include <igl/colon.h>
#include <igl/setdiff.h>
#include <igl/boundary_facets.h>
#include <igl/unique.h>

#include "../include/utils.hpp"


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


typedef Eigen::Triplet<double> T;

void build_laplacian_matrix(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::SparseMatrix<double> &L, std::set<int>& neighbors, int fixed_vertex = -1)
{
    std::vector<std::vector<int>> adj_list;
    igl::adjacency_list(F, adj_list);

    std::vector<T> tripletList;
    
    // On réserve de la place (environ 7 coefficients par sommet en moyenne) -> formule d'euler
    tripletList.reserve(V.rows() * 7);
    
    bool use_neighbors = !neighbors.empty();

    for(int i = 0; i < V.rows(); ++i)
    {

        if (i == fixed_vertex) continue; 

        if (use_neighbors && neighbors.find(i) == neighbors.end()) continue;

        // valence de i (nb voisins)
        double k = (double)adj_list[i].size();

        if(k == 0) {
            tripletList.push_back(T(i, i, 1.0));
            continue;
        }

        tripletList.push_back(T(i, i, 1.0));
        double coeff_k = -1.0 / k;
        for (int voisin : adj_list[i])
        {
            tripletList.push_back(T(i, voisin, coeff_k));
        }
    }
    L.resize(V.rows(), V.rows());
    L.setFromTriplets(tripletList.begin(), tripletList.end());
}

void solve_implicit_heat(const Eigen::VectorXd &U_old, 
                         const Eigen::SparseMatrix<double> &L, 
                         double lambda, 
                         Eigen::VectorXd &U_new)
{
    int n = U_old.rows();

    Eigen::SparseMatrix<double> I(n, n);
    I.setIdentity();

    Eigen::SparseMatrix<double> A = I + lambda * L ;

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

    solver.compute(A);
    if(solver.info() != Eigen::Success) {
        std::cerr << "Erreur de factorisation !" << std::endl;
        return;
    }
    U_new = solver.solve(U_old);
}


void solve_cotangent_heat(const Eigen::VectorXd &U_old, 
                          const Eigen::SparseMatrix<double> &L_cot,
                          const Eigen::SparseMatrix<double> &M,
                          double lambda, 
                          Eigen::VectorXd &U_new)
{
    
    Eigen::SparseMatrix<double> A = M - lambda * L_cot;
    Eigen::VectorXd b = M * U_old;

    // 3. Résolution
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
    
    solver.compute(A);
    if(solver.info() != Eigen::Success) {
        std::cerr << "Erreur de factorisation (Cotangent) !" << std::endl;
        return;
    }

    U_new = solver.solve(b);
}


void solve_laplace_system(const Eigen::SparseMatrix<double> &L_cot,
                        const Eigen::MatrixXd &V, 
                        const Eigen::MatrixXi &F,
                        Eigen::SparseMatrix<double> &L_in_in,
                        Eigen::SparseMatrix<double> &L_in_b,

                        Eigen::VectorXi &all,
                        Eigen::VectorXi &in, 
                        Eigen::VectorXi &b,
                        Eigen::VectorXi &IA,
                        Eigen::VectorXi &IC,
                
                        Eigen::VectorXd &bc,

                        const ViewerControls &ui_state,
                        const std::set<int> &neighbors,
                    
                        Eigen::VectorXd &U_in)
{
    //boundary edges
    std::vector<int> b_vec;
    std::vector<double> bc_vec;

    b_vec.push_back(ui_state.selected_vertex);
    bc_vec.push_back(1.0);

    if (!neighbors.empty()) {
        for(int i=0; i<V.rows(); ++i) {
            if(neighbors.find(i) == neighbors.end()) {
                b_vec.push_back(i);
                bc_vec.push_back(0.0);
            }
        }
    } else {
        if(ui_state.selected_vertex != 0) {
             b_vec.push_back(0); 
             bc_vec.push_back(0.0);
        } else {
             b_vec.push_back(1);
             bc_vec.push_back(0.0);
        }
    }

    b.resize(b_vec.size());
    bc.resize(bc_vec.size());
    for(size_t i=0; i<b_vec.size(); ++i) {
        b(i) = b_vec[i];
        bc(i) = bc_vec[i];
    }

    // List of all vertex indices
    igl::colon<int>(0,V.rows()-1,all);

    // List of interior indices
    igl::setdiff(all,b,in,IA);

    //Laplacian slicing
    igl::slice(L_cot,in,in,L_in_in);
    igl::slice(L_cot,in,b,L_in_b);


    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(-L_in_in);

    if(solver.info() != Eigen::Success) {
        std::cerr << "Erreur de factorisation !" << std::endl;
        return;
    }

    Eigen::VectorXd rhs = L_in_b * bc;
    U_in = solver.solve(rhs);
}