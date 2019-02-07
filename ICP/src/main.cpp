#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <iostream>
//#include "tutorial_shared_path.h"

class MyContext
{
public:

	MyContext():nv_len(0.2), point_size(20),line_width(2){}
	~MyContext() {}

	float nv_len; 
	float point_size;
	float line_width;

	void reset_display(igl::opengl::glfw::Viewer& viewer)
	{
		const Eigen::MatrixXd V = (Eigen::MatrixXd(8, 3) <<
			0.0, 0.0, 0.0,
			0.0, 0.0, 1.0,
			0.0, 1.0, 0.0,
			0.0, 1.0, 1.0,
			1.0, 0.0, 0.0,
			1.0, 0.0, 1.0,
			1.0, 1.0, 0.0,
			1.0, 1.0, 1.0).finished();

		const Eigen::MatrixXi F = (Eigen::MatrixXi(12, 3) <<
			1, 7, 5,
			1, 3, 7,
			1, 4, 3,
			1, 2, 4,
			3, 8, 7,
			3, 4, 8,
			5, 7, 8,
			5, 8, 6,
			1, 5, 6,
			1, 6, 2,
			2, 6, 8,
			2, 8, 4).finished().array() - 1;

		viewer.data().clear();
		viewer.data().set_mesh(V, F);
		viewer.core.align_camera_center(V, F);

		viewer.data().line_width = line_width;
		viewer.data().point_size = point_size;

		viewer.data().add_points(V, Eigen::RowVector3d(1, 1, 0));

		Eigen::MatrixXd V2 = V.array() + nv_len;
		viewer.data().add_edges(V, V2, Eigen::RowVector3d(1, 0, 1));
	}

private:

};

MyContext g_myctx;


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

	// menu variable Shared between two menus
	double doubleVariable = 0.1f; 

	// Add content to the default menu window via defining a Lambda expression with captures by reference([&])
	menu.callback_draw_viewer_menu = [&]()
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
			enum Orientation { Up = 0, Down, Left, Right };
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
			if (num_choices != (int)choices.size())
			{
				choices.resize(num_choices);
				for (int i = 0; i < num_choices; ++i)
					choices[i] = std::string(1, 'A' + i);
				if (idx_choice >= num_choices)
					idx_choice = num_choices - 1;
			}
			ImGui::Combo("Letter", &idx_choice, choices);

			// Add a button
			if (ImGui::Button("Print Hello", ImVec2(-1, 0)))
			{
				std::cout << "Hello\n";
			}
		}
	};

	// Add additional windows via defining a Lambda expression with captures by reference([&])
	menu.callback_draw_custom_window = [&]()
	{
		// Define next window position + size
		ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(300, 160), ImGuiSetCond_FirstUseEver);
		ImGui::Begin( "MyProperties", nullptr, ImGuiWindowFlags_NoSavedSettings );
		
		// point size
		// [event handle] if value changed
		if (ImGui::InputFloat("point_size", &g_myctx.point_size))
		{
			std::cout << "point_size changed\n";
			viewer.data().point_size = g_myctx.point_size;
		}

		// line width
		// [event handle] if value changed
		if(ImGui::InputFloat("line_width", &g_myctx.line_width))
		{
			std::cout << "line_width changed\n";
			viewer.data().line_width = g_myctx.line_width;
		}

		// length of normal line
		// [event handle] if value changed
		if (ImGui::InputFloat("line_length", &g_myctx.nv_len))
		{
			//pass
			std::cout << "line_length changed\n";
			g_myctx.reset_display(viewer);
		}
		
		ImGui::End();
	};

	// registered a event handler
	viewer.callback_key_down = &key_down;

	g_myctx.reset_display(viewer);

	// Call GUI
	viewer.launch();

}
