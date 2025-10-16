#ifndef VIEWER_CONTROLS_H
#define VIEWER_CONTROLS_H

#include <vector>
#include <string>

// Forward declarations pour éviter les dépendances lourdes dans le .h
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
    // Méthode principale de rendu de l'interface
    void render(igl::opengl::glfw::imgui::ImGuiMenu& menu);

private:
    // Énumération pour la direction
    enum Orientation { Up = 0, Down, Left, Right };

    // Variables d'état de l'interface
    double doubleVariable = 0.0;
    bool boolVariable = true;
    Orientation direction = Up;
    int num_choices = 3;
    std::vector<std::string> choices;
    int idx_choice = 0;

    // Méthodes de rendu des widgets individuels
    void render_double_input();
    void render_bool_checkbox();
    void render_direction_combo();
    void render_letter_selector();
    void render_hello_button();

    // Méthode utilitaire
    void update_letter_choices();
};

#endif // VIEWER_CONTROLS_H
