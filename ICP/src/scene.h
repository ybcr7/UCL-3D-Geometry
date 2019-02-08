#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>

class Scene
{
public:
    Scene(igl::opengl::glfw::Viewer& refviewer):viewer(refviewer){}
    ~Scene(){}
    
    // Reset Scene
    void Reset();
    
    // Task 1 & Task 4
    void Point2PointAlign();
    
    // Task 2
    void MatchRotation();
    
    // Task 3
    void AddNoise();
    
    // Task 5
    void MuiltMeshAlign();
    
    // Task 6
    void Point2PlaneAlign();
    
private:
    
    igl::opengl::glfw::Viewer& viewer;
    
    Eigen::MatrixXd V1, V2, V3, V4, V5;
    Eigen::MatrixXi F1, F2, F3, F4, F5;
    
    
    
};
