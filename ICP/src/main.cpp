#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <iostream>
#include "scene.h"

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier){
	std::cout << "Key: " << key << " " << (unsigned int)key << std::endl;
	return false;
}

int main(int argc, char *argv[]){
	// Init the viewer
	igl::opengl::glfw::Viewer viewer;

	// Attach a menu plugin
	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);

	// Menu variable Shared between two menus
	double rotationZValue = 30.0f;
    double gaussianSD = 0.0f;
    int iteration = 50;

    Scene* scene = new Scene(viewer);

    menu.callback_draw_viewer_menu = [&]()
    {
        // Draw parent menu content
        menu.draw_viewer_menu();
        
        // Add a panel
        if (ImGui::CollapsingHeader("Utilities", ImGuiTreeNodeFlags_DefaultOpen))
        {
            if (ImGui::Button("Reset Scene", ImVec2(-1, 0)))
            {
                scene->Reset();
            }
            
            if(ImGui::InputInt("Iteration", &iteration))
            {
                scene->SetIteration(iteration);
            }
        }
    };
    
    menu.callback_draw_custom_window = [&]()
    {
        ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 0), ImGuiSetCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(512, 128), ImGuiSetCond_FirstUseEver);
        ImGui::Begin( "Operations", nullptr, ImGuiWindowFlags_NoSavedSettings );
        
        if (ImGui::CollapsingHeader("Task 1", ImGuiTreeNodeFlags_NoAutoOpenOnLog))
        {
            if (ImGui::Button("Align Meshes", ImVec2(-1, 0))){
                scene->Point2PointAlign();
            }
        }
        
        if (ImGui::CollapsingHeader("Task 2 ~ 4", ImGuiTreeNodeFlags_NoAutoOpenOnLog))
        {
            ImGui::InputDouble("Degree (Z-Axis)", &rotationZValue, 0, 0, "%.4f");
            ImGui::InputDouble("Zero-Mean Gaussian SD", &gaussianSD, 0, 0, "%.4f");
            
            if (ImGui::Button("Rotate Mesh with Noise", ImVec2(-1, 0))){
                scene->RotateMeshWithNoise(rotationZValue,gaussianSD);
            }
            
            if (ImGui::Button("Align Meshes", ImVec2(-1, 0))){
                scene->Point2PointAlign();
            }
            
            if (ImGui::Button("Align Meshes (Optimised)", ImVec2(-1, 0))){
                scene->Point2PointAlignOptimised();
            }
        }
        
        if (ImGui::CollapsingHeader("Task 5", ImGuiTreeNodeFlags_NoAutoOpenOnLog))
        {
            if (ImGui::Button("Align Multiple Meshes", ImVec2(-1, 0))){
                scene->MuiltMeshAlign();
            }
        }
        
        if (ImGui::CollapsingHeader("Task 6", ImGuiTreeNodeFlags_NoAutoOpenOnLog))
        {
            if (ImGui::Button("Align Meshes (Normal-Based)", ImVec2(-1, 0))){
                scene->Point2PlaneAlign();
            }
        }
        
        ImGui::End();
    };

	// Registered a event handler
	viewer.callback_key_down = &key_down;

    // Init the scene
	scene->Reset();

	// Call GUI
	viewer.launch();

}
