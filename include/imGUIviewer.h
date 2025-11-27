#ifndef VIEWER_CONTROLS_H
#define VIEWER_CONTROLS_H

#include <vector>
#include <string>

namespace igl {
    namespace opengl {
        namespace glfw {
            namespace imgui {
                class ImGuiMenu;
            }
        }
    }
}

class ViewerControls {
public:
    void render(igl::opengl::glfw::imgui::ImGuiMenu& menu);
    int selected_vertex = 0;
    int k_rings = 0;
    
    bool boundary = false;
    bool selection_mode = false;

private:
    enum Orientation { Up = 0, Down, Left, Right };

    // Variables
    bool boolVariable = true;
    Orientation direction = Up;
    int num_choices = 3;
    std::vector<std::string> choices;
    int idx_choice = 0;

    void render_k_ring();
    void render_boundary();
    void render_bool_checkbox();
    void render_direction_combo();
    void render_letter_selector();
    void render_hello_button();

    void update_letter_choices();
};

#endif // VIEWER_CONTROLS_H
