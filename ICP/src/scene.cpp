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
    
    Visualise(rendering_data.size());

}

void Scene::Point2PointAlign(){
    
    rendering_data.clear();
    
    Eigen::MatrixXd Vx = V2;
    
    for (size_t i=0; i<iteration;i++){
        
        Vx = ICP::ICPBasic(V1, Vx);
        
        Eigen::MatrixXd V(V1.rows()+Vx.rows(), V1.cols());
        V << V1,Vx;
        
        Eigen::MatrixXi F(F1.rows()+F2.rows(),F1.cols());
        F << F1,(F2.array()+V1.rows());
        Eigen::MatrixXd C(F.rows(),3);
        C <<
        Eigen::RowVector3d(1.0,0.5,0.25).replicate(F1.rows(),1),
        Eigen::RowVector3d(1.0,0.8,0.0).replicate(F2.rows(),1);
        
        if (i+1==iteration){
            
            // Find non-overlapping area
            // Vx to V1
            std::pair<Eigen::MatrixXi, Eigen::MatrixXi> FF2 = ICP::FindNonOverlappingFaces(V1, Vx, F2);
            // V1 to Vx
            std::pair<Eigen::MatrixXi, Eigen::MatrixXi> FF1 = ICP::FindNonOverlappingFaces(Vx, V1, F1);

            // No changes to vertex
            Eigen::MatrixXd V(V1.rows() + V1.rows() + Vx.rows()+ Vx.rows(), V1.cols());
            V << V1,V1,Vx,Vx;
            
            //
            Eigen::MatrixXi F(FF1.first.rows()+FF1.second.rows()+FF2.first.rows()+FF2.second.rows(),F1.cols());
            F << FF1.first, (FF1.second.array()+V1.rows()), (FF2.first.array()+V1.rows()+V1.rows()),(FF2.second.array()+Vx.rows() + V1.rows()+V1.rows());
            
            Eigen::MatrixXd C(F.rows(),3);
            C <<
            Eigen::RowVector3d(1.0,0.8,0.0).replicate(FF1.first.rows(),1),
            Eigen::RowVector3d(1.0,0.0,0.0).replicate(FF1.second.rows(),1),
            Eigen::RowVector3d(1.0,0.5,0.25).replicate(FF2.first.rows(),1),
            Eigen::RowVector3d(1.0,0.0,0.0).replicate(FF2.second.rows(),1);
            
            rendering_data.push_back(RenderingData{V,F,C});
        }else{

            rendering_data.push_back(RenderingData{V,F,C});
        }
        
    }
    
    Visualise(rendering_data.size());

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
    
    Visualise(rendering_data.size());

}

void Scene::Point2PointAlignOptimised(){
    
    
    
    
    
    
    
    
    
}

void Scene::LoadMultiple(){
    
    rendering_data.clear();
    
    igl::readOFF("../data/bun000.off", V1, F1);
    igl::readOFF("../data/bun045.off", V2, F2);
    igl::readOFF("../data/bun315.off", V3, F3);
    igl::readOFF("../data/bun090.off", V4, F4);
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
    
    Visualise(rendering_data.size());
}



void Scene::MuiltMeshAlign(){
    
    Eigen::MatrixXd V1r = V1;
    Eigen::MatrixXi F1r = F1;
    
    Eigen::MatrixXd V2r = V2;
    Eigen::MatrixXi F2r = F2;
    
    Eigen::MatrixXd V3r = V3;
    Eigen::MatrixXi F3r = F3;
    
    Eigen::MatrixXd V4r = V4;
    Eigen::MatrixXi F4r = F4;
    
    Eigen::MatrixXd V5r = V5;
    Eigen::MatrixXi F5r = F5;
    
    for (size_t i=0; i<iteration;i++){
        
        V2r = ICP::ICPBasic(V1, V2r);

        Eigen::MatrixXd Vx(V1r.rows()+V2r.rows(), V1r.cols());
        Vx << V1r,V2r;

        Eigen::MatrixXi Fx(F1r.rows()+F2r.rows(),F1r.cols());
        Fx << F1r, (F2r.array()+V1r.rows());

        Eigen::MatrixXd Cx(Fx.rows(),3);
        Cx <<
        Eigen::RowVector3d(1.0,0.5,0.25).replicate(F1r.rows(),1),
        Eigen::RowVector3d(1.0,0.8,0.0).replicate(F2r.rows(),1);


        V3r = ICP::ICPBasic(Vx, V3r);

        Eigen::MatrixXd Vy(Vx.rows()+V3r.rows(), Vx.cols());
        Vy << Vx,V3r;

        Eigen::MatrixXi Fy(Fx.rows()+F3r.rows(),Fx.cols());
        Fy << Fx, (F3r.array()+Vx.rows());

        Eigen::MatrixXd Cy(Fy.rows(),3);
        Cy <<
        Eigen::RowVector3d(1.0,0.5,0.25).replicate(Fx.rows(),1),
        Eigen::RowVector3d(1.0,0.8,0.0).replicate(F3r.rows(),1);

        
        V4r = ICP::ICPBasic(Vy, V4r);
//
        Eigen::MatrixXd Vz(Vy.rows()+V4r.rows(), Vy.cols());
        Vz << Vy,V4r;
//
        Eigen::MatrixXi Fz(Fy.rows()+F4r.rows(),Fy.cols());
        Fz << Fy, (F4r.array()+Vy.rows());
//
        Eigen::MatrixXd Cz(Fz.rows(),3);
        Cz <<
        Eigen::RowVector3d(1.0,0.5,0.25).replicate(Fy.rows(),1),
        Eigen::RowVector3d(1.0,0.8,0.0).replicate(F4r.rows(),1);

        
        
        //rendering_data.push_back(RenderingData{Vy,Fy,Cy});
        rendering_data.push_back(RenderingData{Vz,Fz,Cz});
        
    }
    
//
//
//        V_to_process = ICP::ICPBasic(V_target, V_to_process);
//
////        Eigen::MatrixXd V(V_target.rows()+V_to_process.rows(), V_target.cols());
////        V <<
////        V_target,
////        V_to_process;
////
////        Eigen::MatrixXi F(F_target.rows()+F_to_process.rows(),F_target.cols());
////        F <<
////        F_target,
////        (F_to_process.array()+V_target.rows());
////
////        Eigen::MatrixXd C(F.rows(),3);
////        C <<
////        Eigen::RowVector3d(1.0,0.5,0.25).replicate(F_target.rows(),1),
////        Eigen::RowVector3d(1.0,0.8,0.0).replicate(F_to_process.rows(),1);
//
//        if (i + 1 == iteration){
//
//            Eigen::MatrixXd V(V_target.rows()+V_to_process.rows(), V_target.cols());
//            V <<
//            V_target,
//            V_to_process;
//
//            Eigen::MatrixXi F(F_target.rows()+F_to_process.rows(),F_target.cols());
//            F <<
//            F_target,
//            (F_to_process.array()+V_target.rows());
//
//            Eigen::MatrixXd C(F.rows(),3);
//            C <<
//            Eigen::RowVector3d(1.0,0.5,0.25).replicate(F_target.rows(),1),
//            Eigen::RowVector3d(1.0,0.8,0.0).replicate(F_to_process.rows(),1);
//
//            rendering_data.push_back(RenderingData{V,F,C});
//
//
//        }
//
//
//    }
//
//    V_target = rendering_data.back().V;
//    F_target = rendering_data.back().F;
//    V_to_process = V3;
//    F_to_process = F3;
//
//    for (size_t i=0; i<iteration;i++){
//
//        V_to_process = ICP::ICPBasic(V_target, V_to_process);
////
////        Eigen::MatrixXd V(V_target.rows()+V_to_process.rows(), V_target.cols());
////        V <<
////        V_target,
////        V_to_process;
////
////        Eigen::MatrixXi F(F_target.rows()+F_to_process.rows(),F_target.cols());
////        F <<
////        F_target,
////        (F_to_process.array()+V_target.rows());
////
////        Eigen::MatrixXd C(F.rows(),3);
////        C <<
////        Eigen::RowVector3d(1.0,0.5,0.25).replicate(F_target.rows(),1),
////        Eigen::RowVector3d(1.0,0.8,0.0).replicate(F_to_process.rows(),1);
//
//
//        if (i+1==iteration){
//
//            Eigen::MatrixXd V(V_target.rows()+V_to_process.rows(), V_target.cols());
//            V <<
//            V_target,
//            V_to_process;
//
//            Eigen::MatrixXi F(F_target.rows()+F_to_process.rows(),F_target.cols());
//            F <<
//            F_target,
//            (F_to_process.array()+V_target.rows());
//
//            Eigen::MatrixXd C(F.rows(),3);
//            C <<
//            Eigen::RowVector3d(1.0,0.5,0.25).replicate(F_target.rows(),1),
//            Eigen::RowVector3d(1.0,0.8,0.0).replicate(F_to_process.rows(),1);
//
//            rendering_data.push_back(RenderingData{V,F,C});
//        }
//
//    }
    
    //cummulated_data = rendering_data.back();

    Visualise(rendering_data.size());
}


void Scene::Point2PlaneAlign(){
    
    
    
}

void Scene::SetIteration(int i){
    if (i < 1){
        iteration = 1;
    }else{
        iteration = i;
    }
}

void Scene::Visualise(int i){
    if (i > 0 && i <= rendering_data.size()){
        viewer.data().clear();
        viewer.data().set_mesh(rendering_data[i-1].V, rendering_data[i-1].F);
        viewer.data().set_colors(rendering_data[i-1].C);
        viewer.data().set_face_based(true);
    }else{
        viewer.data().clear();
    }
}
