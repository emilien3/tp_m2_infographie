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


   // Corners of the bounding box
    // Eigen::MatrixXd V_box(8,3);
    // V_box <<
    // m(0), m(1), m(2),
    // M(0), m(1), m(2),
    // M(0), M(1), m(2),
    // m(0), M(1), m(2),
    // m(0), m(1), M(2),
    // M(0), m(1), M(2),
    // M(0), M(1), M(2),
    // m(0), M(1), M(2);
    
    // // Edges of the bounding box
    // Eigen::MatrixXi E_box(12,2);
    // E_box <<
    // 0, 1,
    // 1, 2,
    // 2, 3,
    // 3, 0,
    // 4, 5,
    // 5, 6,
    // 6, 7,
    // 7, 4,
    // 0, 4,
    // 1, 5,
    // 2, 6,
    // 7 ,3;
    

