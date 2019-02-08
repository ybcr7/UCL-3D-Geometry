#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <iostream>
#include "scene.h"

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
	std::cout << "Key: " << key << " " << (unsigned int)key << std::endl;
	return false;
}

int main(int argc, char *argv[])
{
	// Init the viewer
	igl::opengl::glfw::Viewer viewer;

	// Attach a menu plugin
	igl::opengl::glfw::imgui::ImGuiMenu menu;
    
	viewer.plugins.push_back(&menu);

	// Menu variable Shared between two menus
	double doubleVariable = 0.1f; 

    Scene* scene = new Scene(viewer);
    
	// Add content to the default menu window via defining a Lambda expression with captures by reference([&])
	menu.callback_draw_viewer_menu = [&]()
	{
		// Draw parent menu content
		menu.draw_viewer_menu();
        
        // Add a panel
        if (ImGui::CollapsingHeader("Operations", ImGuiTreeNodeFlags_DefaultOpen))
        {
            if (ImGui::Button("Reset Scene", ImVec2(-1, 0)))
            {
                scene->Reset();
            }
            
            if (ImGui::Button("Point-to-Point", ImVec2(-1, 0)))
            {
                scene->Point2PointAlign();
            }
            
            if (ImGui::Button("Match Rotation", ImVec2(-1, 0)))
            {
                scene->MatchRotation();
            }
            
            if (ImGui::Button("Add Noise", ImVec2(-1, 0)))
            {
                scene->AddNoise();
            }
            
            if (ImGui::Button("Point-to-Point (Optimised)", ImVec2(-1, 0)))
            {
                scene->Point2PointAlignOptimised();
            }
            
            if (ImGui::Button("Multi-Mesh", ImVec2(-1, 0)))
            {
                scene->Reset();
            }
            
            if (ImGui::Button("Point-to-Plane", ImVec2(-1, 0)))
            {
                scene->Point2PlaneAlign();
            }
        }
	};

	// Registered a event handler
	viewer.callback_key_down = &key_down;

    // Init the scene
	scene->Reset();

	// Call GUI
	viewer.launch();

}
