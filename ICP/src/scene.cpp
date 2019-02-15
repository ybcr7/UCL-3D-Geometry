#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include "scene.h"
#include "icp.h"

#define FILE_PATH "../data/"

struct Scene::RenderingData{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXd C;
};

Scene::Scene(igl::opengl::glfw::Viewer& refViewer):viewer(refViewer){
    iteration = 300;
    subsample_rate = 0;
    mark_out = false;
}

Scene::~Scene(){}

void Scene::Initialise(){
    
    rendering_data.clear();

    igl::readOFF(FILE_PATH "bun000.off", V1, F1);
    igl::readOFF(FILE_PATH "bun045.off", V2, F2);

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

        // Basic ICP algorithm
        Vx = ICP::ICPBasic(V1, Vx);

        // Generate data and store them for display
        Eigen::MatrixXd V(V1.rows()+Vx.rows(), V1.cols());
        V << V1,Vx;
        Eigen::MatrixXi F(F1.rows()+F2.rows(),F1.cols());
        F << F1,(F2.array()+V1.rows());
        Eigen::MatrixXd C(F.rows(),3);
        C <<
        Eigen::RowVector3d(1.0,0.5,0.25).replicate(F1.rows(),1),
        Eigen::RowVector3d(1.0,0.8,0.0).replicate(F2.rows(),1);

        // If the V not longer changes or reaches the iteration limit
        if (i+1 == iteration){

            std::cout << "Iteration times: " + std::to_string(i+1) << std::endl;

            if (mark_out){
                // Find non-overlapping area
                // Vx to V1
                std::pair<Eigen::MatrixXi, Eigen::MatrixXi> FF2 = ICP::FindNonOverlappingFaces(V1, Vx, F2);
                // V1 to Vx
                std::pair<Eigen::MatrixXi, Eigen::MatrixXi> FF1 = ICP::FindNonOverlappingFaces(Vx, V1, F1);

                Eigen::MatrixXd V(V1.rows() + V1.rows() + Vx.rows()+ Vx.rows(), V1.cols());
                V << V1,V1,Vx,Vx;

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

            break;

        }else{

            rendering_data.push_back(RenderingData{V,F,C});
        }

    }

    Visualise(rendering_data.size());

}

void Scene::RotateMeshWithNoise(double x, double y, double z, double sd){
    
    rendering_data.clear();
    
    // Load M1
    igl::readOFF(FILE_PATH "bun000.off", V1, F1);
    
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

//    rendering_data.clear();
//
//    Eigen::MatrixXd Vx = V2;
//
//    // Get a subsample from V2
//    Eigen::MatrixXd Vs = ICP::GetSubsample(V2, subsample_rate);
//
//    // Use the subsample to perform ICP algorithm
//    for (size_t i=0; i<iteration;i++){
//
//
//        Eigen::MatrixXd V_matched = ICP::FindCorrespondences(V1, Vs);
//        std::pair<Eigen::Matrix3d, Eigen::RowVector3d> transform = ICP::EstimateRigidTransform(Vs, V_matched);
//        Vs = ICP::ApplyRigidTransform(Vs, transform);
//
//        // Also apply the transform to the original V2 so that we can see how the estimated transforms from subsamples can affect the final result
//        Vx = ICP::ApplyRigidTransform(Vx, transform);
//
//        // Generate data and store them for display
//        Eigen::MatrixXd V(V1.rows()+Vx.rows(), V1.cols());
//        V << V1,Vx;
//        Eigen::MatrixXi F(F1.rows()+F2.rows(),F1.cols());
//        F << F1,(F2.array()+V1.rows());
//        Eigen::MatrixXd C(F.rows(),3);
//        C <<
//          Eigen::RowVector3d(1.0,0.5,0.25).replicate(F1.rows(),1),
//                Eigen::RowVector3d(1.0,0.8,0.0).replicate(F2.rows(),1);
//
//        if (mark_out && i+1==iteration){
//
//            // Find non-overlapping area
//            // Vx to V1
//            std::pair<Eigen::MatrixXi, Eigen::MatrixXi> FF2 = ICP::FindNonOverlappingFaces(V1, Vx, F2);
//            // V1 to Vx
//            std::pair<Eigen::MatrixXi, Eigen::MatrixXi> FF1 = ICP::FindNonOverlappingFaces(Vx, V1, F1);
//
//            Eigen::MatrixXd V(V1.rows() + V1.rows() + Vx.rows()+ Vx.rows(), V1.cols());
//            V << V1,V1,Vx,Vx;
//
//            Eigen::MatrixXi F(FF1.first.rows()+FF1.second.rows()+FF2.first.rows()+FF2.second.rows(),F1.cols());
//            F << FF1.first, (FF1.second.array()+V1.rows()), (FF2.first.array()+V1.rows()+V1.rows()),(FF2.second.array()+Vx.rows() + V1.rows()+V1.rows());
//
//            Eigen::MatrixXd C(F.rows(),3);
//            C <<
//              Eigen::RowVector3d(1.0,0.8,0.0).replicate(FF1.first.rows(),1),
//                    Eigen::RowVector3d(1.0,0.0,0.0).replicate(FF1.second.rows(),1),
//                    Eigen::RowVector3d(1.0,0.5,0.25).replicate(FF2.first.rows(),1),
//                    Eigen::RowVector3d(1.0,0.0,0.0).replicate(FF2.second.rows(),1);
//
//            rendering_data.push_back(RenderingData{V,F,C});
//        }else{
//            rendering_data.push_back(RenderingData{V,F,C});
//        }
//    }
//
//    Visualise(rendering_data.size());

}

void Scene::LoadMultiple(){
    
    rendering_data.clear();
    
    igl::readOFF(FILE_PATH "bun000.off", V1, F1);
    igl::readOFF(FILE_PATH "bun045.off", V2, F2);
    igl::readOFF(FILE_PATH "bun315.off", V3, F3);
    igl::readOFF(FILE_PATH "bun090.off", V4, F4);
    igl::readOFF(FILE_PATH "bun180.off", V5, F5);
    
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

void Scene::MultiMeshAlign(){

    rendering_data.clear();

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

//    for (size_t i=0; i<iteration;i++){
//
//        V2r = ICP::FindBestStartRotation(V1, V2r);
//
//        V2r = ICP::ICPBasic(V1, V2r);
//
//        Eigen::MatrixXd Vx(V1r.rows()+V2r.rows(), V1r.cols());
//        Vx << V1r,V2r;
//
//        Eigen::MatrixXi Fx(F1r.rows()+F2r.rows(),F1r.cols());
//        Fx << F1r, (F2r.array()+V1r.rows());
//
//        Eigen::MatrixXd Cx(Fx.rows(),3);
//        Cx <<
//        Eigen::RowVector3d(1.0,0.5,0.25).replicate(F1r.rows(),1),
//        Eigen::RowVector3d(1.0,0.8,0.0).replicate(F2r.rows(),1);
//
//        V3r = ICP::FindBestStartRotation(Vx, V3r);
//        V3r = ICP::ICPBasic(Vx, V3r);
//
//        Eigen::MatrixXd Vy(Vx.rows()+V3r.rows(), Vx.cols());
//        Vy << Vx,V3r;
//
//        Eigen::MatrixXi Fy(Fx.rows()+F3r.rows(),Fx.cols());
//        Fy << Fx, (F3r.array()+Vx.rows());
//
//        Eigen::MatrixXd Cy(Fy.rows(),3);
//        Cy <<
//        Eigen::RowVector3d(1.0,0.5,0.25).replicate(Fx.rows(),1),
//        Eigen::RowVector3d(1.0,0.8,0.0).replicate(F3r.rows(),1);
//
//        V4r = ICP::FindBestStartRotation(Vy, V4r);
//        V4r = ICP::ICPBasic(Vy, V4r);
//
//        Eigen::MatrixXd Vz(Vy.rows()+V4r.rows(), Vy.cols());
//        Vz << Vy,V4r;
//
//        Eigen::MatrixXi Fz(Fy.rows()+F4r.rows(),Fy.cols());
//        Fz << Fy, (F4r.array()+Vy.rows());
//
//        Eigen::MatrixXd Cz(Fz.rows(),3);
//        Cz <<
//        Eigen::RowVector3d(1.0,0.5,0.25).replicate(Fy.rows(),1),
//        Eigen::RowVector3d(1.0,0.8,0.0).replicate(F4r.rows(),1);
//
//        V5r = ICP::FindBestStartRotation(Vz, V5r);
//        V5r = ICP::ICPBasic(Vz, V5r);
//
//        Eigen::MatrixXd Vo(Vz.rows()+V5r.rows(), Vz.cols());
//        Vo << Vz,V5r;
//
//        Eigen::MatrixXi Fo(Fz.rows()+F5r.rows(),Fz.cols());
//        Fo << Fz, (F5r.array()+Vz.rows());
//
//        Eigen::MatrixXd Co(Fo.rows(),3);
//        Co <<
//        Eigen::RowVector3d(1.0,0.5,0.25).replicate(Fz.rows(),1),
//        Eigen::RowVector3d(1.0,0.8,0.0).replicate(F5r.rows(),1);
//
//
//        if (i+1==iteration){
//            Eigen::MatrixXd C(Fo.rows(),3);
//            C<<
//            Eigen::RowVector3d(1.0,0.5,0.25).replicate(F1r.rows(),1),
//            Eigen::RowVector3d(1.0,0.8,0.0).replicate(F2r.rows(),1),
//            Eigen::RowVector3d(0.25,0.6,1.0).replicate(F3r.rows(),1),
//            Eigen::RowVector3d(0.2,0.7,0.45).replicate(F4r.rows(),1),
//            Eigen::RowVector3d(0.8,0.0,0.8).replicate(F5r.rows(),1);
//            rendering_data.push_back(RenderingData{Vo,Fo,C});
//        }else{
//            rendering_data.push_back(RenderingData{Vo,Fo,Co});
//        }
//
//    }



//
//    V2r = ICP::FindBestStartRotation(V1, V2r);
//    for (size_t i=0; i<iteration;i++) {
//
//        Eigen::MatrixXd V_matched = ICP::FindCorrespondences(V1, V2r);
//        std::pair<Eigen::Matrix3d, Eigen::RowVector3d> transform = ICP::EstimateRigidTransform(V_matched,V2r);
//        V2r = ICP::ApplyRigidTransform(V2r, transform);
//
//        if (transform.second.norm() < loop_threshold || i + 1 == iteration){
//            std::cout << "Iteration times: " + std::to_string(i+1) << std::endl;
//            break;
//        }
//    }
//
//    Eigen::MatrixXd V12(V1.rows()+V2r.rows(), V1.cols());
//    V12<<V1, V2r;
//    Eigen::MatrixXi F12(F1.rows()+F2r.rows(), F1.cols());
//    F12<<F1, (F2r.array()+V1.rows());
//
//    V3r = ICP::FindBestStartRotation(V2r, V3r);
//    for (size_t i=0; i<iteration;i++) {
//
//        Eigen::MatrixXd V_matched = ICP::FindCorrespondences(V12, V3r);
//        std::pair<Eigen::Matrix3d, Eigen::RowVector3d> transform = ICP::EstimateRigidTransform(V_matched,V3r);
//        V3r = ICP::ApplyRigidTransform(V3r, transform);
//
//        if (transform.second.norm() < loop_threshold || i + 1 == iteration){
//            std::cout << "Iteration times: " + std::to_string(i+1) << std::endl;
//            break;
//        }
//    }
//
//    Eigen::MatrixXd V123(V12.rows()+V3r.rows(), V1.cols());
//    V123<<V12, V3r;
//    Eigen::MatrixXi F123(F12.rows()+F3r.rows(), F1.cols());
//    F123<<F12, (F3r.array()+V12.rows());
//
//    V4r = ICP::FindBestStartRotation(V3r, V4r);
//    for (size_t i=0; i<iteration;i++) {
//
//        Eigen::MatrixXd V_matched = ICP::FindCorrespondences(V123, V4r);
//        std::pair<Eigen::Matrix3d, Eigen::RowVector3d> transform = ICP::EstimateRigidTransform(V_matched,V4r);
//        V4r = ICP::ApplyRigidTransform(V4r, transform);
//
//        if (transform.second.norm() < loop_threshold || i + 1 == iteration){
//            std::cout << "Iteration times: " + std::to_string(i+1) << std::endl;
//            break;
//        }
//    }
//
//    Eigen::MatrixXd V1234(V123.rows()+V4r.rows(), V1.cols());
//    V1234<<V123, V4r;
//    Eigen::MatrixXi F1234(F123.rows()+F4r.rows(), F1.cols());
//    F1234<<F123,(F4r.array()+V123.rows());
//
//    V5r = ICP::FindBestStartRotation(V4r, V5r);
//    for (size_t i=0; i<iteration;i++) {
//
//        Eigen::MatrixXd V_matched = ICP::FindCorrespondences(V1234, V5r);
//        std::pair<Eigen::Matrix3d, Eigen::RowVector3d> transform = ICP::EstimateRigidTransform(V_matched,V5r);
//        V5r = ICP::ApplyRigidTransform(V5r, transform);
//
//        if (transform.second.norm() < loop_threshold || i + 1 == iteration){
//            std::cout << "Iteration times: " + std::to_string(i+1) << std::endl;
//            break;
//        }
//    }
//
//    Eigen::MatrixXd V12345(V1234.rows()+V5.rows(), V1.cols());
//    V12345<<V1234, V5r;
//    Eigen::MatrixXi F12345(F1234.rows()+F5r.rows(), F1.cols());
//    F12345<<F1234,(F5r.array()+V1234.rows());
//
//    Eigen::MatrixXd C(F12345.rows(),3);
//    C<<
//    Eigen::RowVector3d(1.0,0.5,0.25).replicate(F1.rows(),1),
//    Eigen::RowVector3d(1.0,0.8,0.0).replicate(F2r.rows(),1),
//    Eigen::RowVector3d(0.25,0.6,1.0).replicate(F3r.rows(),1),
//    Eigen::RowVector3d(0.2,0.7,0.45).replicate(F4r.rows(),1),
//    Eigen::RowVector3d(0.8,0.0,0.8).replicate(F5r.rows(),1);
//
//    rendering_data.push_back(RenderingData{V12345,F12345,C});
//
//    Visualise(rendering_data.size());
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

void Scene::SetMarkOut(bool b) {
    mark_out = b;
}

void Scene::Visualise(int i){
    if (i > 0 && i <= rendering_data.size()){
        viewer.data().clear();
        viewer.data().set_mesh(rendering_data[i-1].V, rendering_data[i-1].F);
        viewer.data().set_colors(rendering_data[i-1].C);
        viewer.data().set_face_based(true);
    }else {
        viewer.data().clear();
    }
}

void Scene::SetSubsampleRate(double s){
    if (s >= 0 && s <= 99){
        subsample_rate = s;
    }else {
        subsample_rate = 0;
    }
    std::cout << "Subsample Rate: " + std::to_string(subsample_rate) + "%" << std::endl;
}
