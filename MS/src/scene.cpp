#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/file_exists.h>
#include "scene.h"
#include "ms.h"

Scene::Scene(igl::opengl::glfw::Viewer& refViewer):viewer(refViewer){
    default_C << 1.0,1.0,0.0;
}

Scene::~Scene(){}

void Scene::Discretisation(int mode){
    Eigen::VectorXd C_out;
    switch (mode){
        case 0:
            C_out = MS::UniformMeanCurvature(V,F);

            break;
        case 1:
            C_out = MS::GaussianCurvature(V,F);

            break;
        case 2:
            C_out = MS::NonUniformMeanCurvature(V,F);

            break;
        default:
            std::cout << "ERROR: Undefined Discretisation Mode" << std::endl;
            break;
    }

    C_out = 1000 * C_out.array() / (C_out.maxCoeff() - C_out.minCoeff());
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
            V_out = MS::ExplicitSmoothing(V_out,F,lambda,iteration);
            C_out = MS::UniformMeanCurvature(V_out,F);
            C_out = 5 * C_out.array() / (C_out.maxCoeff() - C_out.minCoeff());
            igl::parula(C_out, false, C);
            Visualise(V_out,F);
            break;
        case 1:
            V_out = MS::ImplicitSmoothing(V_out,F,lambda,iteration);
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

    std::vector<const char *> FILE_PATH_LIST{
            "../data/",
            "../../data/",
            "../../../data/"};
    bool file_found = false;

    for (const auto & FILE_PATH : FILE_PATH_LIST)
    {
        if ( igl::file_exists(FILE_PATH + filename) )
        {
            if (igl::readOFF(FILE_PATH+filename, V, F)) {
                file_found = true;
                break;
            }
        }
    }

    if (!file_found) {
        std::cout << "ERROR: "<< filename << " not found" << std::endl;
        exit(1);
    }

    C.resize(V.rows(),V.cols());
    for (int i = 0; i < V.rows(); i++){
        C.row(i) = default_C;
    }
    Visualise(V,F);
}

void Scene::AddNoise(){
    V = MS::AddNoise(V, noise);
    Visualise(V, F);
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

void Scene::Visualise(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in){
    viewer.data().clear();
    viewer.data().show_overlay_depth = 1;
    viewer.data().set_mesh(V_in, F_in);
    viewer.data().set_colors(C);
    viewer.core.align_camera_center(V_in, F_in);
}
