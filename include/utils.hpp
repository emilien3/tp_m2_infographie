#include <Eigen/Sparse>
#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <sstream>
#include <set>
#include <vector>
#include "../include/imGUIviewer.h"


void boundingBox(Eigen::MatrixXd &V_box);

std::set<int>& get_k_ring_neighbors(std::set<int>& neighbors, 
                const std::vector<std::vector<int>>& adjacency_list,
                int vertex_id, 
                int k);

void anneaux_topologiques(std::set<int>& neighbors, Eigen::MatrixXd& V,
                          Eigen::MatrixXi& F,
                          Eigen::VectorXd& neighborhood_colors, 
                          int selected_vertex, 
                          int k_rings);

bool callback_visualize(igl::opengl::glfw::Viewer& viewer, 
                        Eigen::MatrixXd& C, 
                        Eigen::MatrixXd& V, 
                        Eigen::MatrixXi& F,
                        const Eigen::VectorXd& values);

bool callback_colorize_face(igl::opengl::glfw::Viewer& viewer, 
                            Eigen::MatrixXd& C, 
                            Eigen::MatrixXd& V, 
                            Eigen::MatrixXi& F);

bool callback_select_mode(igl::opengl::glfw::Viewer& viewer, 
                            Eigen::MatrixXd& V,
                            Eigen::MatrixXi& F, 
                            ViewerControls& ui_state);

void compute_mean_curvature(Eigen::VectorXd& s_laplacien, 
                            Eigen::MatrixXd& V, 
                            Eigen::MatrixXi& F);

void step_heat_diffusion(Eigen::VectorXd& U, 
                        const std::vector<std::vector<int>>& adj_list, 
                        std::set<int>& neighbors);

void build_laplacian_matrix(const Eigen::MatrixXd &V, 
                            const Eigen::MatrixXi &F, 
                            Eigen::SparseMatrix<double> &L, 
                            std::set<int>& neighbors,
                            int fixed_vertex);

void solve_implicit_heat(const Eigen::VectorXd &U_old, 
                         const Eigen::SparseMatrix<double> &L, 
                         double lambda, 
                         Eigen::VectorXd &U_new);

void solve_cotangent_heat(const Eigen::VectorXd &U_old, 
                            const Eigen::SparseMatrix<double> &L_cot,
                            const Eigen::SparseMatrix<double> &M,
                            double lambda, 
                            Eigen::VectorXd &U_new);


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
                        const std::set<int>& neighbors,
                        
                        Eigen::VectorXd &U_in);
