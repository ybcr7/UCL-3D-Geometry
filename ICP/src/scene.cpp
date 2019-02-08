#include "scene.h"

void Scene::Reset(){
    
    viewer.data().clear();

    igl::readOFF("../data/bun000.off", V1, F1);
    igl::readOFF("../data/bun180.off", V2, F2);

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

void Scene::MatchRotation(){
    
    
    
}


void Scene::AddNoise(){
    
    
    
    
}


void Scene::MuiltMeshAlign(){
//    viewer.data().clear;
//
//    igl::readOFF(TUTORIAL_SHARED_PATH "/bun000.off", V1, F1);
//    igl::readOFF(TUTORIAL_SHARED_PATH "/bun045.off", V2, F2);
//    igl::readOFF(TUTORIAL_SHARED_PATH "/bun090.off", V3, F3);
//    igl::readOFF(TUTORIAL_SHARED_PATH "/bun180.off", V4, F4);
//    igl::readOFF(TUTORIAL_SHARED_PATH "/bun270.off", V5, F5);
    
}


void Scene::Point2PlaneAlign(){
    
    
    
    
    
    
    
}


