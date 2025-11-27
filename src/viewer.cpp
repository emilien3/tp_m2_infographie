#include "../include/imGUIviewer.h"
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <imgui.h>
#include <iostream>
#include <algorithm>

void ViewerControls::render(igl::opengl::glfw::imgui::ImGuiMenu& menu) {
    menu.draw_viewer_menu();
    
    if (ImGui::CollapsingHeader("Laplacien config", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Checkbox("Selection Mode", &selection_mode);

        render_k_ring();
        render_boundary();
        render_bool_checkbox();
        render_direction_combo();
        render_letter_selector();
        render_hello_button();
    }
}

void ViewerControls::render_k_ring() {
    ImGui::InputInt("k_rings", &k_rings);
}


void ViewerControls::render_boundary() {
    ImGui::Checkbox("boundary", &boundary);
}

void ViewerControls::render_bool_checkbox() {
    // if (ImGui::Checkbox("bool", &boolVariable)) {
    //     std::cout << "boolVariable: " << std::boolalpha << boolVariable << std::endl;
    // }
}

void ViewerControls::render_direction_combo() {
    // ImGui::Combo("Direction", (int*)(&direction), "Up\0Down\0Left\0Right\0\0");
}

void ViewerControls::render_letter_selector() {
    // if (ImGui::InputInt("Num letters", &num_choices)) {
    //     num_choices = std::max(1, std::min(26, num_choices));
    // }
    
    // if (num_choices != (int)choices.size()) {
    //     update_letter_choices();
    // }
    
    // ImGui::Combo("Letter", &idx_choice, choices);
}

void ViewerControls::update_letter_choices() {
    // choices.resize(num_choices);
    // for (int i = 0; i < num_choices; ++i) {
    //     choices[i] = std::string(1, 'A' + i);
    // }
    // if (idx_choice >= num_choices) {
    //     idx_choice = num_choices - 1;
    // }
}

void ViewerControls::render_hello_button() {
    // if (ImGui::Button("Print Hello", ImVec2(-1, 0))) {
    //     std::cout << "Hello\n";
    // }
}
