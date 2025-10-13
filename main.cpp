#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <iostream>
#include <sstream>

Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd C;


int main(int argc, char *argv[])
{

    // viewer.load_mesh_from_file(".../data/sphere.obj");

    // Load a mesh in OFF format
    igl::readOFF("../data/bunny.off", V, F);
    // Save the mesh in OBJ format
    // igl::writeOBJ("cube.obj",V,F);

    Eigen::Vector3d m = V.colwise().minCoeff();
    Eigen::Vector3d M = V.colwise().maxCoeff();

    // Corners of the bounding box
    Eigen::MatrixXd V_box(8,3);
    V_box <<
    m(0), m(1), m(2),
    M(0), m(1), m(2),
    M(0), M(1), m(2),
    m(0), M(1), m(2),
    m(0), m(1), M(2),
    M(0), m(1), M(2),
    M(0), M(1), M(2),
    m(0), M(1), M(2);
    // Edges of the bounding box
    Eigen::MatrixXi E_box(12,2);
    E_box <<
    0, 1,
    1, 2,
    2, 3,
    3, 0,
    4, 5,
    5, 6,
    6, 7,
    7, 4,
    0, 4,
    1, 5,
    2, 6,
    7 ,3;


    // Plot the mesh //
    igl::opengl::glfw::Viewer viewer;
    // Parametrage du viewer 
    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
    // Setting du viewer
    // V for vertices and F for faces

    //// PLOT THE BOUNDING BOX 
    // Plot the corners of the bounding box as points
    viewer.data().add_points(V_box,Eigen::RowVector3d(1,0,0));

    // Plot the edges of the bounding box
    for (unsigned i=0;i<E_box.rows(); ++i)
        viewer.data().add_edges
        (
        V_box.row(E_box(i,0)),
        V_box.row(E_box(i,1)),
        Eigen::RowVector3d(1,0,0)
        );

    // Plot labels with the coordinates of bounding box vertices
    std::stringstream l1;
    l1 << m(0) << ", " << m(1) << ", " << m(2);
    viewer.data().add_label(m+Eigen::Vector3d(-0.007, 0, 0),l1.str());
    std::stringstream l2;
    l2 << M(0) << ", " << M(1) << ", " << M(2);
    viewer.data().add_label(M+Eigen::Vector3d(0.007, 0, 0),l2.str());
    // activate label rendering
    viewer.data().show_custom_labels = true;

    /////////////////
    ///// COLOR /////
    ////////////////


    // Set the colors for the model based on the cartesien colors

     // Use the (normalized) vertex positions as colors
    C =
        (V.rowwise()            - V.colwise().minCoeff()).array().rowwise()/
        (V.colwise().maxCoeff() - V.colwise().minCoeff()).array();

    // Add per-vertex colors
    viewer.data().set_colors(C);





    /////////////////
    ///// ImGUI /////
    ////////////////

    double doubleVariable = 0.1f;

    // Rendering of text labels is handled by ImGui, so we need to enable the ImGui
    // plugin to show text labels.
    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);
    menu.callback_draw_viewer_window = [&]()
    {
        // Draw parent menu content
        menu.draw_viewer_menu();

        // Add new group
        if (ImGui::CollapsingHeader("New Group", ImGuiTreeNodeFlags_DefaultOpen))
        {
        // Expose variable directly ...
        ImGui::InputDouble("double", &doubleVariable, 0, 0, "%.4f");

        // ... or using a custom callback
        static bool boolVariable = true;
        if (ImGui::Checkbox("bool", &boolVariable))
        {
            // do something
            std::cout << "boolVariable: " << std::boolalpha << boolVariable << std::endl;
        }

        // Expose an enumeration type
        enum Orientation { Up=0, Down, Left, Right };
        static Orientation dir = Up;
        ImGui::Combo("Direction", (int *)(&dir), "Up\0Down\0Left\0Right\0\0");

        // We can also use a std::vector<std::string> defined dynamically
        static int num_choices = 3;
        static std::vector<std::string> choices;
        static int idx_choice = 0;
        if (ImGui::InputInt("Num letters", &num_choices))
        {
            num_choices = std::max(1, std::min(26, num_choices));
        }
        if (num_choices != (int) choices.size())
        {
            choices.resize(num_choices);
            for (int i = 0; i < num_choices; ++i)
            choices[i] = std::string(1, 'A' + i);
            if (idx_choice >= num_choices)
            idx_choice = num_choices - 1;
        }
        ImGui::Combo("Letter", &idx_choice, choices);

        // Add a button
        if (ImGui::Button("Print Hello", ImVec2(-1,0)))
        {
            std::cout << "Hello\n";
        }
        }
    };


    // Draw additional windows
    menu.callback_draw_custom_window = [&]()
    {
        // Define next window position + size
        ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiCond_FirstUseEver);
        ImGui::Begin(
            "New Window", nullptr,
            ImGuiWindowFlags_NoSavedSettings
        );

        // Expose the same variable directly ...
        ImGui::PushItemWidth(-80);
        ImGui::DragScalar("double", ImGuiDataType_Double, &doubleVariable, 0.1, 0, 0, "%.4f");
        ImGui::PopItemWidth();

        static std::string str = "bunny";
        ImGui::InputText("Name", str);

        ImGui::End();
    };


    viewer.data().set_mesh(V, F);

    viewer.data().add_label(viewer.data().V.row(0) + viewer.data().V_normals.row(0).normalized()*0.005, "Hello World!");

    viewer.launch();
    /////////////////////////
        

}
