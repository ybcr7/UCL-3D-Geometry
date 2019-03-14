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
    std::pair<Eigen::MatrixXd, Eigen::MatrixXi> mesh_out;
    switch (mode){
        case 0:
            mesh_out = MS::UniformMeanCurvature(V,F);
            break;
        case 1:
            mesh_out = MS::UniformGaussianCurvature(V,F);
            break;
        case 2:
            mesh_out = MS::NonUniformMeanCurvature(V,F);
            break;
        default:
            std::cout << "ERROR: Undefined Discretisation Mode" << std::endl;
            break;
    }
}


void Scene::Reconstruction() {

}


void Scene::Smoothing(int mode) {
    std::pair<Eigen::MatrixXd, Eigen::MatrixXi> mesh_out;
    switch (mode){
        case 0:
            mesh_out = MS::ExplicitSmoothing(V,F);
            break;
        case 1:
            mesh_out = MS::ImplicitSmoothing(V,F);
            break;
        default:
            std::cout << "ERROR: Undefined Smoothing Mode" << std::endl;
            break;
    }
}

void Scene::Initialise(std::string filename){

    rendering_data.clear();

    igl::readOFF(FILE_PATH + filename, V, F);

    Eigen::MatrixXd Vx(V.rows(), V.cols());
    Vx << V;
    Eigen::MatrixXi Fx(F.rows(),F.cols());
    Fx << F;
    Eigen::MatrixXd Cx(F.rows(),3);
    Cx <<
    Eigen::RowVector3d(1.0,0.5,0.25).replicate(F.rows(),1),

    rendering_data.push_back(RenderingData{Vx,Fx,Cx});

    Visualise(rendering_data.size());

}





void Scene::AddNoise(){

//    rendering_data.clear();
//
//    // Load M1
//    igl::readOFF(FILE_PATH "bun000.off", V1, F1);
//    igl::readOFF(FILE_PATH "bun045.off", V2, F2);
//
//    // M2' = M2
//    V2 = ICP::AddNoise(V2, sd);
//
//    // Display meshes
//    Eigen::MatrixXd V(V1.rows()+V2.rows(), V1.cols());
//    V << V1,V2;
//
//    Eigen::MatrixXi F(F1.rows()+F2.rows(),F1.cols());
//    F << F1, (F2.array()+V1.rows());
//    Eigen::MatrixXd C(F.rows(),3);
//    C <<
//      Eigen::RowVector3d(1.0,0.5,0.25).replicate(F1.rows(),1),
//            Eigen::RowVector3d(1.0,0.8,0.0).replicate(F2.rows(),1);
//
//    rendering_data.push_back(RenderingData{V,F,C});
//
//    Visualise(rendering_data.size());

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



void Scene::Visualise(int i){
    if (i > 0 && i <= rendering_data.size()){
        viewer.data().clear();
        viewer.data().set_mesh(rendering_data[i-1].V, rendering_data[i-1].F);
        viewer.data().set_colors(rendering_data[i-1].C);
        viewer.data().set_face_based(true);
        viewer.core.align_camera_center(rendering_data[i-1].V, rendering_data[i-1].F);
    }else {
        viewer.data().clear();
    }
}
