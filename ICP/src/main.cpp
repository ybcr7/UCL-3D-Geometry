#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <iostream>
#include "scene.h"

int main(int argc, char *argv[]){
	
    // Initialise the viewer
	igl::opengl::glfw::Viewer viewer;

	// Attach a menu plugin
	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);
    
    Scene scene(viewer);

	// Menu variable Shared between two menus
	double rotation_x = 0.0;
    double rotation_y = 0.0;
    double rotation_z = 0.0;
    double gaussian_sd = 0.0;
    int iteration = 50;

    // Draw an optional panel for adjusting global variables
    menu.callback_draw_viewer_menu = [&]()
    {
        // Draw parent menu content
        menu.draw_viewer_menu();
        
        // Add a panel
        if (ImGui::CollapsingHeader("Utilities", ImGuiTreeNodeFlags_DefaultOpen))
        {
            if (ImGui::Button("Reset Scene", ImVec2(-1, 0)))
            {
                scene.Initialise();
            }
            
            if(ImGui::InputInt("Iteration", &iteration))
            {
                scene.SetIteration(iteration);
            }
        }
    };
    
    // Draw a separate panel for performing tasks
    menu.callback_draw_custom_window = [&]()
    {
        ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 0), ImGuiSetCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(512, 128), ImGuiSetCond_FirstUseEver);
        ImGui::Begin( "Operations", nullptr, ImGuiWindowFlags_NoSavedSettings );
        
        if (ImGui::CollapsingHeader("Task 1", ImGuiTreeNodeFlags_NoAutoOpenOnLog))
        {
            if (ImGui::Button("Align Meshes", ImVec2(-1, 0))){
                scene.Point2PointAlign();
            }
        }
        
        if (ImGui::CollapsingHeader("Task 2 ~ 4", ImGuiTreeNodeFlags_NoAutoOpenOnLog))
        {
            ImGui::InputDouble("Degree (X-Axis)", &rotation_x, 0, 0, "%.4f");
            ImGui::InputDouble("Degree (Y-Axis)", &rotation_y, 0, 0, "%.4f");
            ImGui::InputDouble("Degree (Z-Axis)", &rotation_z, 0, 0, "%.4f");
            ImGui::InputDouble("Zero-Mean Gaussian SD", &gaussian_sd, 0, 0, "%.4f");
            
            if (ImGui::Button("Rotate Mesh with Noise", ImVec2(-1, 0))){
                scene.RotateMeshWithNoise(rotation_x,rotation_y,rotation_z,gaussian_sd);
            }
            
            if (ImGui::Button("Align Meshes", ImVec2(-1, 0))){
                scene.Point2PointAlign();
            }
            
            if (ImGui::Button("Align Meshes (Optimised)", ImVec2(-1, 0))){
                scene.Point2PointAlignOptimised();
            }
        }
        
        if (ImGui::CollapsingHeader("Task 5", ImGuiTreeNodeFlags_NoAutoOpenOnLog))
        {
            if (ImGui::Button("Align Multiple Meshes", ImVec2(-1, 0))){
                scene.MuiltMeshAlign();
            }
        }
        
        if (ImGui::CollapsingHeader("Task 6", ImGuiTreeNodeFlags_NoAutoOpenOnLog))
        {
            if (ImGui::Button("Align Meshes (Normal-Based)", ImVec2(-1, 0))){
                scene.Point2PlaneAlign();
            }
        }
        
        ImGui::End();
    };

    // Initialise the scene
	scene.Initialise();

	// Call GUI
	viewer.launch();

}
