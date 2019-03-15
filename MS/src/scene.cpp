#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include "scene.h"
#include "ms.h"

#define FILE_PATH "data/"

struct Scene::RenderingData{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXd C;
};

Scene::Scene(igl::opengl::glfw::Viewer& refViewer):viewer(refViewer){
    iteration = 100;
    lambda = 0;
    noise = 0;
}

Scene::~Scene(){}

void Scene::Discretisation(int mode){
    Eigen::VectorXd C_out;
    switch (mode){
        case 0:
            C_out = MS::UniformMeanCurvature(V,F);
            break;
        case 1:
            C_out = MS::UniformGaussianCurvature(V,F);
            break;
        case 2:
            C_out = MS::NonUniformMeanCurvature(V,F);
            break;
        default:
            std::cout << "ERROR: Undefined Discretisation Mode" << std::endl;
            break;
    }

    C_out = 5 * C_out.array() / (C_out.maxCoeff() - C_out.minCoeff());
    igl::parula(C_out, false, C);

    Visualise();
}

void Scene::Reconstruction() {
    MS::BarycentricArea(V,F);
}


void Scene::Smoothing(int mode) {
    std::pair<Eigen::MatrixXd, Eigen::MatrixXi> mesh_out;
    switch (mode){
        case 0:
            //mesh_out = MS::ExplicitSmoothing(V,F);
            break;
        case 1:
            //mesh_out = MS::ImplicitSmoothing(V,F);
            break;
        default:
            std::cout << "ERROR: Undefined Smoothing Mode" << std::endl;
            break;
    }
}

void Scene::Initialise(std::string filename){

    igl::readOFF(FILE_PATH + filename, V, F);

    Visualise();
}

void Scene::AddNoise(){

    rendering_data.clear();

    V = MS::AddNoise(V, noise);

    Visualise();

}

void Scene::SetNumEigenvector(int e) {
    if (e < 1) {
        iteration = 1;
    } else {
        iteration = e;
    }
}

void Scene::SetIteration(int i) {
    if (i < 1) {
        iteration = 1;
    } else {
        iteration = i;
    }
}

void Scene::SetLambda(double l) {
    if (l < 0) {
        lambda = 0;
    } else {
        lambda = l;
    }
}

void Scene::SetNoise(double n) {
    if (n < 0) {
        noise = 0;
    } else {
        noise = n;
    }
}

void Scene::Visualise(){
    viewer.data().clear();
    viewer.data().show_lines = 0;
    viewer.data().show_overlay_depth = 1;
    viewer.data().set_mesh(V, F);
    viewer.data().set_colors(C);
    viewer.core.align_camera_center(V, F);
}
