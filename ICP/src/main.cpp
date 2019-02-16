#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "scene.h"

int main(int argc, char *argv[]){

    srand (time(NULL));

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
    int iteration = 300;
    double subsample_rate = 0.0;
    int frame = 1;
    bool mark_out = false;

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

            if (ImGui::Checkbox("Show Non-Overlapping Area", &mark_out))
            {

                scene.SetMarkOut(mark_out);

                if (mark_out){

                    scene.Visualise(2);

                }else{

                    scene.Visualise(1);

                }
            }

//            if(ImGui::InputInt("Frame", &frame))
//            {
//                scene.Visualise(frame);
//            }
        }
    };
    
    // Draw a separate panel for performing tasks
    menu.callback_draw_custom_window = [&]()
    {
        ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 0), ImGuiSetCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(448, 384), ImGuiSetCond_FirstUseEver);
        ImGui::Begin( "Operations", nullptr, ImGuiWindowFlags_NoSavedSettings );

        if (ImGui::CollapsingHeader("Task 1", ImGuiTreeNodeFlags_DefaultOpen))
        {
            if (ImGui::Button("Align Meshes", ImVec2(-1, 0))){
                scene.Point2PointAlign();
            }
        }
        
        if (ImGui::CollapsingHeader("Task 2", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::InputDouble("Degree (X-Axis)", &rotation_x, 0, 0, "%.4f");
            ImGui::InputDouble("Degree (Y-Axis)", &rotation_y, 0, 0, "%.4f");
            ImGui::InputDouble("Degree (Z-Axis)", &rotation_z, 0, 0, "%.4f");

            if (ImGui::Button("Rotate Mesh", ImVec2(-1, 0))){
                scene.RotateMesh(rotation_x,rotation_y,rotation_z);
            }
            
            if (ImGui::Button("Align Meshes", ImVec2(-1, 0))){
                scene.Point2PointAlign();
            }
        }

        if (ImGui::CollapsingHeader("Task 3", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::InputDouble("Zero-Mean Gaussian SD", &gaussian_sd, 0, 0, "%.4f");

            if (ImGui::Button("Add Noise", ImVec2(-1, 0))){
                scene.AddNoiseToMesh(gaussian_sd);
            }

            if (ImGui::Button("Align Meshes", ImVec2(-1, 0))){
                scene.Point2PointAlign();
            }

        }

        if (ImGui::CollapsingHeader("Task 4", ImGuiTreeNodeFlags_DefaultOpen))
        {

            if(ImGui::InputDouble("Subsample Rate %", &subsample_rate, 0, 0, "%.4f"))
            {
                scene.SetSubsampleRate(subsample_rate);
            }

            if (ImGui::Button("Align Meshes (Optimised)", ImVec2(-1, 0))){
                scene.Point2PointAlignOptimised();
            }
        }
        
        if (ImGui::CollapsingHeader("Task 5", ImGuiTreeNodeFlags_DefaultOpen))
        {
            if (ImGui::Button("Load Multiple Meshes", ImVec2(-1, 0))){
                scene.LoadMultiple();
            }
            
            if (ImGui::Button("Align Multiple Meshes", ImVec2(-1, 0))){
                scene.MultiMeshAlign();
            }
        }
        
        if (ImGui::CollapsingHeader("Task 6", ImGuiTreeNodeFlags_DefaultOpen))
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