#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include "scene.h"
#include "icp.h"

Scene::Scene(igl::opengl::glfw::Viewer& refViewer):viewer(refViewer){
    iteration = 50;
}

Scene::~Scene(){}

void Scene::Reset(){
    
    viewer.data().clear();

    igl::readOFF("../data/bun000.off", V1, F1);
    igl::readOFF("../data/bun045.off", V2, F2);

    Eigen::MatrixXd V(V1.rows()+V2.rows(), V1.cols());
    V << V1,V2;

    Eigen::MatrixXi F(F1.rows()+F2.rows(),F1.cols());
    F << F1, (F2.array()+V1.rows());
    Eigen::MatrixXd C(F.rows(),3);
    C <<
    Eigen::RowVector3d(0.66,0.0,0.0).replicate(F1.rows(),1),
    Eigen::RowVector3d(0.0,0.0,0.66).replicate(F2.rows(),1);

    viewer.data().set_mesh(V, F);
    viewer.data().set_colors(C);
    viewer.data().set_face_based(true);
    viewer.core.align_camera_center(V, F);

}

void Scene::Point2PointAlign(){
    
    
    
    
}

void Scene::RotateMeshWithNoise(double degreeZ, double sd){
    
    viewer.data().clear();
    
    // Load M1
    igl::readOFF("../data/bun000.off", V1, F1);
    
    // M2 = R(M1), vertex positions are changed while the face relation remains
    V2 = ICP::Rotate(V1, degreeZ);
    F2 = F1;

    // M2' = M2
    V2 = ICP::AddNoise(V2, sd);
    
    // Display meshes
    Eigen::MatrixXd V(V1.rows()+V2.rows(), V1.cols());
    V << V1,V2;

    Eigen::MatrixXi F(F1.rows()+F2.rows(),F1.cols());
    F << F1, (F2.array()+V1.rows());
    Eigen::MatrixXd C(F.rows(),3);
    C <<
    Eigen::RowVector3d(0.66,0.0,0.0).replicate(F1.rows(),1),
    Eigen::RowVector3d(0.0,0.0,0.66).replicate(F2.rows(),1);

    viewer.data().set_mesh(V, F);
    viewer.data().set_colors(C);
    viewer.data().set_face_based(true);
    viewer.core.align_camera_center(V, F);

}

void Scene::Point2PointAlignOptimised(){
    
    
    
    
    
}


void Scene::MuiltMeshAlign(){
    
    viewer.data().clear();

    igl::readOFF("../data/bun000.off", V1, F1);
    igl::readOFF("../data/bun045.off", V2, F2);
    igl::readOFF("../data/bun090.off", V3, F3);
    igl::readOFF("../data/bun180.off", V4, F4);
    igl::readOFF("../data/bun270.off", V5, F5);
    
    viewer.data().set_mesh(V1, F1);
    viewer.core.align_camera_center(V1, F1);
}


void Scene::Point2PlaneAlign(){
    
    
    
    
    
    
    
}



void Scene::DisplayMeshes(){
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    
//    for (int i = 0; i < Vs.size(); i++ ){
//
//
//
//
//
//
//
//
//
//
//
//    }
//
//
//
//    Eigen::MatrixXd V(V1.rows()+V2.rows(), V1.cols());
//    V << V1,V2;
//
//    Eigen::MatrixXi F(F1.rows()+F2.rows(),F1.cols());
//    F << F1, (F2.array()+V1.rows());
//    Eigen::MatrixXd C(F.rows(),3);
//    C <<
//    Eigen::RowVector3d(0.66,0.0,0.0).replicate(F1.rows(),1),
//    Eigen::RowVector3d(0.0,0.0,0.66).replicate(F2.rows(),1);
//
//    viewer.data().set_mesh(V, F);
//    viewer.data().set_colors(C);
//    viewer.data().set_face_based(true);
//    viewer.core.align_camera_center(V, F);
    
}

void Scene::SetIteration(int i){
    if (i < 1){
        iteration = 1;
    }else{
        iteration = i;
    }
    std::cout<< iteration << std::endl;
}
