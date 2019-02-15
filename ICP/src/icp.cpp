#include <Eigen/Dense>
#include <Eigen/SVD>
#include <random>
#include <iostream>
#include <math.h>
#include <time.h>
#include <igl/boundary_loop.h>
#include <igl/bounding_box.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/fit_plane.h>
#include "nanoflann.hpp"
#include "icp.h"

Eigen::MatrixXd ICP::GetSubsample(Eigen::MatrixXd V_to_process, double subsample_rate){

    // A. Random Sampling


    // B. Uniform Subsampling
    std::vector<int> V_index;
    std::random_device rand;
    std::mt19937 generate(rand());
    std::uniform_int_distribution<> distribution(0, V_to_process.rows()-1);
    int vertex_to_keep = V_to_process.rows() - round(V_to_process.rows()*(subsample_rate/100));

    for (int i = 0; i < vertex_to_keep; i++){
        V_index.push_back(distribution(generate));
    }

    std::sort(V_index.begin(),V_index.end());
    V_index.erase(std::unique(V_index.begin(), V_index.end()), V_index.end());

    Eigen::MatrixXd V_out(V_index.size(), 3);
    V_out.setZero();

    for (int i = 0; i < V_index.size(); i ++){
        V_out.row(i) = V_to_process.row(V_index[i]);
    }

    std::cout << V_out.rows() << std::endl;

    return V_out;
}

Eigen::MatrixXd ICP::GetVertexNormal(Eigen::MatrixXd V_target, Eigen::MatrixXd V_to_process){


    // p <-matched q<-to process
    Eigen::MatrixXd VN_out;
    VN_out.resize(V_to_process.rows(),V_to_process.cols());
    VN_out.setZero();

    const size_t num_result = 8; // 4 or 8
    const size_t max_leaf = 10;

    Eigen::RowVector3d center;
    center = V_to_process.colwise().sum()/ double(V_to_process.rows());


    // Generate the KD tree with nanoflann
    nanoflann::KDTreeEigenMatrixAdaptor<Eigen::MatrixXd> kd_tree_index(V_to_process, max_leaf);
    kd_tree_index.index->buildIndex();

    // For each vertex
    for (size_t v=0; v<VN_out.rows(); v++) {

        // Pick the current vertex for query
        Eigen::RowVector3d query_vertex = V_to_process.row(v);

        // Create a query object
        std::vector<size_t> indexes(num_result);
        std::vector<double> dists_sqr(num_result);

        // Find the closest 1 vertex
        nanoflann::KNNResultSet<double> result(num_result);
        result.init(indexes.data(), dists_sqr.data());
        kd_tree_index.index->findNeighbors(result, query_vertex.data(), nanoflann::SearchParams(max_leaf));

        // Assign founded vertex to output matrix
        Eigen::MatrixXd V_founded(num_result, 3);
        for (size_t i = 0; i < num_result; i++){
            V_founded.row(i) = V_to_process.row(indexes[i]);
        }

        Eigen::RowVector3d N, C;
        igl::fit_plane(V_founded, N, C);

        VN_out.row(v) = N.row(v);

        // Flip normal if pointing to the wrong direction
        if((center(0,0)-V_to_process(v,0)) * VN_out(v,0) + (center(0,1)-V_to_process(v,1)) * VN_out(v,1) + (center(0,2)-V_to_process(v,2)) * VN_out(v,2) > 0) {
            VN_out.row(v) = -VN_out.row(v);
        }
    }

    return VN_out;

}

std::pair<Eigen::MatrixXi, Eigen::MatrixXi> ICP::FindNonOverlappingFaces(Eigen::MatrixXd V_target, Eigen::MatrixXd V_to_process, Eigen::MatrixXi F_to_process){
    
    // Initialise output matrix
    std::vector<int> V_distant;
    std::vector<Eigen::Vector3i> F_raw;
    Eigen::MatrixXi F_non_overlap(0,3);
    Eigen::MatrixXi F_overlap(0,3);
    Eigen::Vector3i empty (-1,-1,-1);
    
    const size_t num_result = 1;
    const size_t max_leaf = 10;
    const double threshold = 0.00001;

    // Generate the KD tree with nanoflann
    nanoflann::KDTreeEigenMatrixAdaptor<Eigen::MatrixXd> kd_tree_index(V_target, max_leaf);
    kd_tree_index.index->buildIndex();

    // For each vertex
    for (size_t v=0; v<V_to_process.rows(); v++){

        // Pick the current vertex for query
        Eigen::RowVector3d query_vertex = V_to_process.row(v);

        // Create a query object
        std::vector<size_t> indexes(num_result);
        std::vector<double> dists_sqr(num_result);

        // Find the closest 1 vertex
        nanoflann::KNNResultSet<double> result(num_result);
        result.init(indexes.data(), dists_sqr.data());
        kd_tree_index.index->findNeighbors(result, query_vertex.data(), nanoflann::SearchParams(max_leaf));

        // Check if it is distant
        if (dists_sqr[0] > threshold){
            V_distant.push_back(v);
        }
    }

    for (size_t f=0; f<F_to_process.rows(); f++){
        F_raw.push_back(F_to_process.row(f));
    }

    // Construct the non-overlapping face list
    for (size_t f=0; f<F_raw.size(); f++){
        
        for (size_t i=0; i < V_distant.size(); i ++){
            
            for (size_t k=0; k<F_raw[f].size(); k++){
                
                if (F_raw[f][k] == V_distant[i]){
                    
                    F_non_overlap.conservativeResize(F_non_overlap.rows()+1, 3);
                    
                    F_non_overlap.row(F_non_overlap.rows()-1) = F_to_process.row(f);
                    
                    F_raw[f] = empty;
                    
                    break;
                }
            }
        }
    }

    // Construct the overlapping face list
    for (size_t i=0; i <F_raw.size(); i++){
        if (F_raw[i] != empty){
            F_overlap.conservativeResize(F_overlap.rows()+1, 3);
            F_overlap.row(F_overlap.rows()-1) = F_raw[i];
        }
    }

    // Output
    return std::pair<Eigen::MatrixXi, Eigen::MatrixXi>(F_overlap, F_non_overlap);
    
}

Eigen::MatrixXd ICP::Rotate(Eigen::MatrixXd V_in, double x, double y, double z){
    
    // Initialise
    Eigen::MatrixXd V_out;
    V_out.resize(V_in.rows(), V_in.cols());
    V_out.setZero();
    
    // Construct rotation matrix
    Eigen::Matrix3d R;
    R = Eigen::AngleAxisd(x * M_PI/180, Eigen::Vector3d::UnitX()) * Eigen::AngleAxisd(y * M_PI/180, Eigen::Vector3d::UnitY()) * Eigen::AngleAxisd(z * M_PI/180, Eigen::Vector3d::UnitZ());
    
    // Rotate the vertex
    V_out = V_in * R;
    
    return V_out;
    
}

Eigen::MatrixXd ICP::AddNoise(Eigen::MatrixXd V_in, double sd){
    
    // Initialise output matrix
    Eigen::MatrixXd V_out;
    V_out.resize(V_in.rows(), V_in.cols());
    V_out.setZero();

    Eigen::MatrixXd V_bounding;
    Eigen::MatrixXi F_bounding;

    igl::bounding_box(V_in, V_bounding, F_bounding);
    double noise_scale_x = abs(V_bounding.row(0).x()-V_bounding.row(4).x());
    double noise_scale_y = abs(V_bounding.row(0).y()-V_bounding.row(2).y());
    double noise_scale_z = abs(V_bounding.row(0).z()-V_bounding.row(1).z());

    // Random number generation
    std::default_random_engine rnd;
    std::normal_distribution<double> gaussian(0.0, sd);
    
    // Add noise to the vertex
    for (int i=0; i<V_out.rows(); i++){
        double x =gaussian(rnd)/(10000*noise_scale_x);
        double y =gaussian(rnd)/(10000*noise_scale_y);
        double z =gaussian(rnd)/(10000*noise_scale_z);
        Eigen::RowVector3d noise(x,y,z);
        V_out.row(i) = V_in.row(i) + noise;
    }
    
    return V_out;
    
}

Eigen::MatrixXd ICP::FindBestStartRotation(Eigen::MatrixXd V_target, Eigen::MatrixXd V_to_process){
    const int rotate_degree = 120;
    std::vector<Eigen::MatrixXd> V_rotate_list;
    std::vector<double> distance_list;
    //Eigen::RowVector3d center_target = V_target.colwise().sum()/V_target.rows();
    // Apply a 0(360), 120 , 240 rotation for each axis (not at the same time)
    for (int x = 0; x < 3; x++){
        for (int y = 0; y < 3; y++){
            for (int z = 0; z < 3; z ++){
                Eigen::MatrixXd V_rotated = Rotate(V_to_process, x*120, y*120, z*120);
                V_rotate_list.push_back(V_rotated);
                Eigen::MatrixXd V_matched = FindCorrespondences(V_target, V_rotated).first;
                Eigen::RowVector3d center_matched = V_matched.colwise().sum()/V_matched.rows();
                Eigen::RowVector3d center_rotated = V_rotated.colwise().sum()/V_rotated.rows();
                double distance = (center_matched-center_rotated).norm();
                distance_list.push_back(distance);
            }
        }
    }

    // Find the one with the rotated one with smallest euclidean distance
    size_t min = std::min_element(distance_list.begin(),distance_list.end()) - distance_list.begin();
    return V_rotate_list[min];
}

//std::pair<Eigen::MatrixXd, Eigen::MatrixXd> ICP::RejectPairs(Eigen::MatrixXd V_matched, Eigen::MatrixXd V_to_process){
//
//    const double k = 2.0;
//
//    Eigen::MatrixXd V_matched_rejected(0,3);
//    Eigen::MatrixXd V_raw_reduced(0,3);
//
//    std::vector<int> error_list;
//
//    Eigen::RowVector3d mean_distance_matched = V_matched.colwise().sum()/V_matched.rows();
//    Eigen::RowVector3d mean_distance_raw = V_to_process.colwise().sum()/V_to_process.rows();
//    double distance = (mean_distance_matched-mean_distance_raw).norm();
//
//    for (size_t i = 0; i < V_matched.rows(); i ++){
//        if ((V_matched.row(i) - V_to_process.row(i)).norm() > k * distance){
//            error_list.push_back(i);
//        }
//    }
//
//    for (size_t i = 0; i < V_matched.rows(); i ++){
//        if(std::find(error_list.begin(), error_list.end(), i) != error_list.end()) {
//            // If i is in the error list
//        } else {
//            // Otherwise assign non-error vertex to outputs
//            V_matched_rejected.conservativeResize(V_matched_rejected.rows()+1, 3);
//            V_matched_rejected.row(V_matched_rejected.rows()-1) = V_matched.row(i);
//
//            V_raw_reduced.conservativeResize(V_raw_reduced.rows()+1, 3);
//            V_raw_reduced.row(V_raw_reduced.rows()-1) = V_to_process.row(i);
//        }
//    }
//
//    std::cout << V_matched_rejected.rows() << std::endl;
//
//    return std::pair<Eigen::MatrixXd, Eigen::MatrixXd>(V_matched_rejected, V_raw_reduced);
//}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> ICP::FindCorrespondences(Eigen::MatrixXd V_target, Eigen::MatrixXd V_to_process){

    // Initialise output matrix
    Eigen::MatrixXd V_out;
    std::vector<double> distances, distances_raw;
    V_out.resize(V_to_process.rows(), V_to_process.cols());
    V_out.setZero();

    std::vector<int> refined_index;

    const double k = 2.0;
    const size_t num_result = 1;
    const size_t max_leaf = 20;

    double distance_median;

    // Generate the KD tree with nanoflann
    nanoflann::KDTreeEigenMatrixAdaptor<Eigen::MatrixXd> kd_tree_index(V_target, max_leaf);
    kd_tree_index.index->buildIndex();

    // For each vertex
    for (size_t v=0; v<V_out.rows(); v++){

        // Pick the current vertex for query
        Eigen::RowVector3d query_vertex = V_to_process.row(v);

        // Create a query object
        std::vector<size_t> indexes(num_result);
        std::vector<double> dists_sqr(num_result);

        // Find the closest 1 vertex
        nanoflann::KNNResultSet<double> result(num_result);
        result.init(indexes.data(), dists_sqr.data());
        kd_tree_index.index->findNeighbors(result, query_vertex.data(), nanoflann::SearchParams(max_leaf));

        // Assign founded vertex to output matrix
        V_out.row(v) = V_target.row(indexes[0]);
        distances.push_back(dists_sqr[0]);

    }

    distances_raw = distances;

    // The pairs rejection can be done using k * median distance
    // Reference slide: http://resources.mpi-inf.mpg.de/deformableShapeMatching/EG2011_Tutorial/slides/2.1%20Rigid%20ICP.pdf

    // Find median distance value
    std::sort(distances.begin(),distances.end());

    if (distances.size() % 2 == 0){
        distance_median = (distances[distances.size()/2-1]+distances[distances.size()/2])/2;
    }else{
        distance_median = distances[distances.size()/2];
    }

    for (int i = 0; i < distances_raw.size(); i++){
        if (distances_raw[i] < k * distance_median){
            refined_index.push_back(i);
        }else{
            // Distant vertex, ignore.
        }
    }

    Eigen::MatrixXd V_refined_out(refined_index.size(), 3);
    Eigen::MatrixXd V_refined_raw(refined_index.size(), 3);

    for (int i = 0; i < refined_index.size(); i++){
        V_refined_out.row(i) = V_out.row(refined_index[i]);
        V_refined_raw.row(i) = V_to_process.row(refined_index[i]);
    }

    std::cout << V_refined_out.rows() << std::endl;

    return std::pair<Eigen::MatrixXd, Eigen::MatrixXd>(V_refined_out, V_refined_raw);
}


std::pair<Eigen::Matrix3d, Eigen::RowVector3d> ICP::EstimateRigidTransform(Eigen::MatrixXd V_matched, Eigen::MatrixXd V_to_process){

    // The rigid transform can be estimated from min(R,t) Sigma_i ||p_i - R*q_i - t||^2
    // t = p_bar - R*q_bar
    // R can be estimated from min(R) Sigma_i ||p_hat_i - R* q_hat_i||^2

    std::pair<Eigen::Matrix3d, Eigen::RowVector3d> transform;

    // Define Barycenters
    // p_bar = 1/m sigma p_i => average
    Eigen::RowVector3d p_bar = V_matched.colwise().mean();
    Eigen::RowVector3d q_bar = V_to_process.colwise().mean();

    Eigen::MatrixXd p_hat = V_matched.rowwise() - p_bar;
    Eigen::MatrixXd q_hat = V_to_process.rowwise() - q_bar;

    // Construct A
    Eigen::Matrix3d A;
    A.setZero();

    for (size_t i=0; i<V_matched.rows(); i++){

        Eigen::Vector3d q_i = q_hat.row(i);
        Eigen::Vector3d p_i = p_hat.row(i);

        A += q_i * p_i.transpose();
    }

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd R = svd.matrixU() * svd.matrixV().transpose();
    Eigen::RowVector3d T = p_bar - q_bar * R;

    transform.first = R;
    transform.second = T;

    return transform;

}

Eigen::MatrixXd ICP::ApplyRigidTransform(Eigen::MatrixXd V_to_process, std::pair<Eigen::Matrix3d, Eigen::RowVector3d> transform){

    Eigen::MatrixXd V_out;
    V_out.resize(V_to_process.rows(), V_to_process.cols());
    V_out.setZero();

    // According to the formula, p = Rq + t
    for (size_t i=0;i<V_to_process.rows();i++){
        Eigen::Vector3d row = V_to_process.row(i);
        V_out.row(i) = (row.transpose() * transform.first + transform.second).transpose();
    }

    return V_out;
}

Eigen::MatrixXd ICP::ICPBasic(Eigen::MatrixXd V_target, Eigen::MatrixXd V_to_process){
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> correspondences = FindCorrespondences(V_target, V_to_process);
    std::pair<Eigen::Matrix3d, Eigen::RowVector3d> transform = ICP::EstimateRigidTransform(correspondences.first, correspondences.second);
    return ICP::ApplyRigidTransform(V_to_process, transform);
}


Eigen::MatrixXd ICP::ICPOptimised(Eigen::MatrixXd V_target, Eigen::MatrixXd V_to_process, double subsample_rate){
    Eigen::MatrixXd V_subsampled = GetSubsample(V_to_process, subsample_rate);
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> correspondences = FindCorrespondences(V_target, V_subsampled);
    std::pair<Eigen::Matrix3d, Eigen::RowVector3d> transform = ICP::EstimateRigidTransform(correspondences.first, correspondences.second);
    return ICP::ApplyRigidTransform(V_to_process, transform);
}

//Eigen::MatrixXd ICP::ICPAdvanced()
