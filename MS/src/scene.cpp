#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include "scene.h"
#include "ms.h"

#define FILE_PATH "data/"

Scene::Scene(igl::opengl::glfw::Viewer& refViewer):viewer(refViewer){

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

    C_out = 100 * C_out.array() / (C_out.maxCoeff() - C_out.minCoeff());
    igl::parula(C_out, false, C);

    Visualise(V,F);
}

void Scene::Reconstruction() {
    Eigen::MatrixXd V_out(V.rows(),V.cols());
    V_out = MS::Reconstruction(V,F,eigenvector);
    Visualise(V_out, F);
}


void Scene::Smoothing(int mode) {
    Eigen::MatrixXd V_out = V;
    Eigen::VectorXd C_out;
    switch (mode){
        case 0:
            for (int i = 0; i < iteration; i ++){
                V_out = MS::ExplicitSmoothing(V_out,F,lambda);
                std::cout << "Complete iteration:" << i << std::endl;
            }
            C_out = MS::UniformMeanCurvature(V_out,F);
            C_out = 5 * C_out.array() / (C_out.maxCoeff() - C_out.minCoeff());
            igl::parula(C_out, false, C);
            Visualise(V_out,F);
            break;
        case 1:
            for (int i = 0; i < iteration; i ++){
                V_out = MS::ImplicitSmoothing(V_out,F,lambda);
                std::cout << "Complete iteration:" << i << std::endl;
            }
            C_out = MS::UniformMeanCurvature(V_out,F);
            C_out = 5 * C_out.array() / (C_out.maxCoeff() - C_out.minCoeff());
            igl::parula(C_out, false, C);
            Visualise(V_out,F);
            break;
        default:
            std::cout << "ERROR: Undefined Smoothing Mode" << std::endl;
            break;
    }
}

void Scene::Initialise(std::string filename){

    igl::readOFF(FILE_PATH + filename, V, F);
    //C = Eigen::Vector3d(1,0,1.0,0.0);
    Visualise(V,F);
}

void Scene::AddNoise(){
    Eigen::MatrixXd V_out = MS::AddNoise(V, noise);
    Visualise(V_out, F);
}

void Scene::SetNumEigenvector(int e) {
    if (e < 1) {
        eigenvector = 1;
    } else {
        eigenvector = e;
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

//void Scene::Visualise(){
//    viewer.data().clear();
//    viewer.data().show_lines = 0;
//    viewer.data().show_overlay_depth = 1;
//    viewer.data().set_mesh(V, F);
//    viewer.data().set_colors(C);
//    viewer.core.align_camera_center(V, F);
//}

void Scene::Visualise(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in){
    viewer.data().clear();
    viewer.data().show_lines = 0;
    viewer.data().show_overlay_depth = 1;
    viewer.data().set_mesh(V_in, F_in);
    viewer.data().set_colors(C);
    viewer.core.align_camera_center(V_in, F_in);
}