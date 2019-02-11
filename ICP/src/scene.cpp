#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include "scene.h"
#include "icp.h"

struct Scene::RenderingData{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXd C;
};

Scene::Scene(igl::opengl::glfw::Viewer& refViewer):viewer(refViewer){
    iteration = 30;
}

Scene::~Scene(){}

void Scene::Initialise(){
    
    rendering_data.clear();

    igl::readOFF("../data/bun000.off", V1, F1);
    igl::readOFF("../data/bun045.off", V2, F2);

    Eigen::MatrixXd V(V1.rows()+V2.rows(), V1.cols());
    V << V1,V2;

    Eigen::MatrixXi F(F1.rows()+F2.rows(),F1.cols());
    F << F1,(F2.array()+V1.rows());
    Eigen::MatrixXd C(F.rows(),3);
    C <<
    Eigen::RowVector3d(1.0,0.5,0.25).replicate(F1.rows(),1),
    Eigen::RowVector3d(1.0,0.8,0.0).replicate(F2.rows(),1);

    rendering_data.push_back(RenderingData{V,F,C});
    
    VisualiseData(rendering_data.size());

}

void Scene::Point2PointAlign(){
    
    rendering_data.clear();
    
    Eigen::MatrixXd Vx = V2;
    
    for (size_t i=0; i<iteration;i++){
        
        Vx = ICP::ICPBasic(V1, Vx, iteration);
        
        Eigen::MatrixXd V(V1.rows()+Vx.rows(), V1.cols());
        V << V1,Vx;
        
        Eigen::MatrixXi F(F1.rows()+F2.rows(),F1.cols());
        F << F1,(F2.array()+V1.rows());
        Eigen::MatrixXd C(F.rows(),3);
        C <<
        Eigen::RowVector3d(1.0,0.5,0.25).replicate(F1.rows(),1),
        Eigen::RowVector3d(1.0,0.8,0.0).replicate(F2.rows(),1);
        
        rendering_data.push_back(RenderingData{V,F,C});
        
    }
    
    VisualiseData(rendering_data.size());

}

void Scene::RotateMeshWithNoise(double x, double y, double z, double sd){
    
    rendering_data.clear();
    
    // Load M1
    igl::readOFF("../data/bun000.off", V1, F1);
    
    // M2 = R(M1), vertex positions are changed while the face relation remains
    V2 = ICP::Rotate(V1, x, y, z);
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
    Eigen::RowVector3d(1.0,0.5,0.25).replicate(F1.rows(),1),
    Eigen::RowVector3d(1.0,0.8,0.0).replicate(F2.rows(),1);

    rendering_data.push_back(RenderingData{V,F,C});
    
    VisualiseData(rendering_data.size());

}

void Scene::Point2PointAlignOptimised(){
    
    

    
}


void Scene::MuiltMeshAlign(){
    
    rendering_data.clear();

    igl::readOFF("../data/bun000.off", V1, F1);
    igl::readOFF("../data/bun045.off", V2, F2);
    igl::readOFF("../data/bun090.off", V3, F3);
    igl::readOFF("../data/bun180.off", V4, F4);
    igl::readOFF("../data/bun270.off", V5, F5);
    
    // Display meshes
    Eigen::MatrixXd V(V1.rows()+V2.rows()+V3.rows()+V4.rows()+V5.rows(), V1.cols());
    V << V1,V2,V3,V4,V5;
    
    Eigen::MatrixXi F(F1.rows()+F2.rows()+F3.rows()+F4.rows()+F5.rows(), F1.cols());
    F <<F1,(F2.array()+V1.rows()), (F3.array()+V2.rows()+V1.rows()), (F4.array()+V3.rows()+V2.rows()+V1.rows()), (F5.array()+V4.rows()+V3.rows()+V2.rows()+V1.rows());
    
    Eigen::MatrixXd C(F.rows(),3);
    C <<
    Eigen::RowVector3d(1.0,0.5,0.25).replicate(F1.rows(),1),
    Eigen::RowVector3d(1.0,0.8,0.0).replicate(F2.rows(),1),
    Eigen::RowVector3d(0.25,0.6,1.0).replicate(F3.rows(),1),
    Eigen::RowVector3d(0.2,0.7,0.45).replicate(F4.rows(),1),
    Eigen::RowVector3d(0.8,0.0,0.8).replicate(F5.rows(),1);
    
    rendering_data.push_back(RenderingData{V,F,C});
    
    VisualiseData(rendering_data.size());
}


void Scene::Point2PlaneAlign(){
    
    
    
    
    
    
    
}

void Scene::SetIteration(int i){
    if (i < 1){
        iteration = 1;
    }else{
        iteration = i;
    }
    std::cout<< iteration << std::endl;
}

void Scene::VisualiseData(int i){
    if (i > 0 && i <= rendering_data.size()){
        viewer.data().clear();
        viewer.data().set_mesh(rendering_data[i-1].V, rendering_data[i-1].F);
        viewer.data().set_colors(rendering_data[i-1].C);
        viewer.data().set_face_based(true);
        viewer.core.align_camera_center(rendering_data[i-1].V, rendering_data[i-1].F);
    }
}
