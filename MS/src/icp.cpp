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

    // #Input: V2, Int
    // #Output: V2_Subsampled

    // Random
    std::vector<int> V_index;
    for (size_t i = 0; i < V_to_process.rows(); i++) {
        if (rand() / double(RAND_MAX) >= subsample_rate/100) {
            V_index.push_back(i);
        }
    }

    Eigen::MatrixXd V_out(V_index.size(), 3);
    V_out.setZero();

    for (size_t i = 0; i < V_index.size(); i++) {
        V_out.row(i) = V_to_process.row(V_index[i]);
    }

    return V_out;

    // Uniform
//    std::vector<int> V_index;
//    std::random_device rand;
//    std::mt19937 generate(rand());
//    std::uniform_int_distribution<> distribution(0, V_to_process.rows() - 1);
//    int vertex_to_keep = round(V_to_process.rows() * (1 - (subsample_rate / 100)));
//
//    for (int i = 0; i < vertex_to_keep; i++) {
//        V_index.push_back(distribution(generate));
//    }
//
//    std::sort(V_index.begin(), V_index.end());
//    V_index.erase(std::unique(V_index.begin(), V_index.end()), V_index.end());
//
//    Eigen::MatrixXd V_out(V_index.size(), 3);
//    V_out.setZero();
//
//    for (int i = 0; i < V_index.size(); i++) {
//        V_out.row(i) = V_to_process.row(V_index[i]);
//    }
//
//    return V_out;

}

Eigen::MatrixXd ICP::GetVertexNormal(Eigen::MatrixXd V_target){

    // #Input: V1
    // #Output: N1

    // p <-matched q<-to process
    Eigen::MatrixXd VN_out;
    VN_out.resize(V_target.rows(),V_target.cols());
    VN_out.setZero();

    const size_t num_result = 20; // 4 or 8
    const size_t max_leaf = 50;

    //Eigen::RowVector3d center;
    Eigen::MatrixXd center(1,3);
    center = V_target.colwise().sum()/ double(V_target.rows());

    // Generate the KD tree with nanoflann
    nanoflann::KDTreeEigenMatrixAdaptor<Eigen::MatrixXd> kd_tree_index(V_target, max_leaf);
    kd_tree_index.index->buildIndex();

    // For each vertex
    for (size_t v=0; v<VN_out.rows(); v++) {

        // Pick the current vertex for query
        Eigen::RowVector3d query_vertex = V_target.row(v);

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
            V_founded.row(i) = V_target.row(indexes[i]);
        }

        Eigen::RowVector3d N, C;
        igl::fit_plane(V_founded, N, C);

        VN_out.row(v) = N;

        // Flip normal if pointing to the wrong direction
        if((center(0,0)-V_target(v,0)) * VN_out(v,0) + (center(0,1)-V_target(v,1)) * VN_out(v,1) + (center(0,2)-V_target(v,2)) * VN_out(v,2) > 0) {
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

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> ICP::FindCorrespondences(Eigen::MatrixXd V_target, Eigen::MatrixXd V_to_process){

    // #Input: V1, V2 (Without Rejection)
    // #Output: V1_Matched, V2_Matched (With Rejection)

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

    //std::cout << "Unprocessed size:" + std::to_string(V_to_process.rows()) << std::endl;


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
    // Reference slide: http://resources.mpi-inf.mpg.de/deformableShapeMatching/EG2011_Tutorial/slides/2.1%20Rigid%20ICP.pdf Page 8

    // Find median distance value
    std::sort(distances.begin(),distances.end());

    if (distances.size() % 2 == 0){
        distance_median = (distances[distances.size()/2-1]+distances[distances.size()/2])/2;
    }else{
        distance_median = distances[distances.size()/2];
    }

    for (int i = 0; i < distances_raw.size(); i++){
        if (distances_raw[i] <= k * distance_median){
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

    //std::cout << "V_refined" + std::to_string(V_refined_out.rows()) << std::endl;

    return std::pair<Eigen::MatrixXd, Eigen::MatrixXd>(V_refined_out, V_refined_raw);

    //return std::pair<Eigen::MatrixXd, Eigen::MatrixXd>(V_out, V_to_process);
}

std::pair<std::pair<Eigen::MatrixXd, Eigen::MatrixXd>, Eigen::MatrixXd> ICP::FindCorrespondencesNormalBased(Eigen::MatrixXd V_target, Eigen::MatrixXd V_to_process, Eigen::MatrixXd N_target){

    // #Input: V1, V2, N1 (Without Rejection)
    // #Output: V1_Matched, N1_Matched, V2_Matched (With Rejection)

    // Initialise output matrix
    Eigen::MatrixXd V_out, N_out;
    std::vector<double> distances, distances_raw;
    V_out.resize(V_to_process.rows(), V_to_process.cols());
    V_out.setZero();

    N_out.resize(V_to_process.rows(), V_to_process.cols());
    N_out.setZero();

    std::vector<int> refined_index;

    const double k = 1.0;
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
    N_out.row(v) = N_target.row(indexes[0]);
    distances.push_back(dists_sqr[0]);

    }

    distances_raw = distances;

    // The pairs rejection can be done using k * median distance
    // Reference slide: http://resources.mpi-inf.mpg.de/deformableShapeMatching/EG2011_Tutorial/slides/2.1%20Rigid%20ICP.pdf Page 8

    // Find median distance value
    std::sort(distances.begin(),distances.end());

    if (distances.size() % 2 == 0){
    distance_median = (distances[distances.size()/2-1]+distances[distances.size()/2])/2;
    }else{
    distance_median = distances[distances.size()/2];
    }

    for (int i = 0; i < distances_raw.size(); i++){
    if (distances_raw[i] <= k * distance_median){
    refined_index.push_back(i);
    }else{
    // Distant vertex, ignore.
    }
    }

    Eigen::MatrixXd V_refined_out(refined_index.size(), 3);
    Eigen::MatrixXd V_refined_raw(refined_index.size(), 3);
    Eigen::MatrixXd N_refined_out(refined_index.size(), 3);

    for (int i = 0; i < refined_index.size(); i++){
    V_refined_out.row(i) = V_out.row(refined_index[i]);
    V_refined_raw.row(i) = V_to_process.row(refined_index[i]);
    N_refined_out.row(i) = N_out.row(refined_index[i]);
    }

    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> refined_target (V_refined_out, N_refined_out);

    return std::pair<std::pair<Eigen::MatrixXd, Eigen::MatrixXd>, Eigen::MatrixXd> (refined_target, V_refined_raw);

}

std::pair<double, std::pair<Eigen::Matrix3d, Eigen::RowVector3d>> ICP::EstimateRigidTransform(Eigen::MatrixXd V_matched, Eigen::MatrixXd V_to_process){

    // #Input: V1_Matched, V2_Matched
    // #Output: R, T

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

        Eigen::Vector3d p_i = p_hat.row(i);
        Eigen::Vector3d q_i = q_hat.row(i);

        A += q_i * p_i.transpose();
    }

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd R = svd.matrixV() * svd.matrixU().transpose();
    Eigen::RowVector3d T = p_bar - (R * q_bar.transpose()).transpose();

    transform.first = R;
    transform.second = T;

    //
    double error_metric = GetErrorMetric(V_matched, ApplyRigidTransform(V_to_process,transform));

    return std::pair<double, std::pair<Eigen::Matrix3d, Eigen::RowVector3d>> (error_metric, transform);

}

std::pair<Eigen::Matrix3d, Eigen::RowVector3d> ICP::EstimateRigidTransformNormalBased(Eigen::MatrixXd V_matched, Eigen::MatrixXd V_to_process, Eigen::MatrixXd N_to_process){

    // #Input: V1_Matched, V2_Matched, N1_Matched
    // #Output: R, T

    // Reference slide: http://resources.mpi-inf.mpg.de/deformableShapeMatching/EG2011_Tutorial/slides/2.1%20Rigid%20ICP.pdf Page 12

    // Construct A and b to solve the point-to-plane error metric
    // A = p_i x n_i, n_i
    // b = -(p_i - q_i) . n_i

    Eigen::MatrixXd A (V_matched.rows(), 6);
    Eigen::MatrixXd b (V_matched.rows(), 1);

    for (size_t i = 0; i < V_matched.rows(); i++){

        Eigen::RowVector3d N = N_to_process.row(i);
        Eigen::RowVector3d S = V_to_process.row(i);
        Eigen::RowVector3d D = V_matched.row(i);

        A(i,0) = N.z()*D.y()-N.y()*D.z();
        A(i,1) = N.x()*D.z()-N.z()*D.x();
        A(i,2) = N.y()*D.x()-N.x()*D.y();
        A(i,3) = N.x();
        A(i,4) = N.y();
        A(i,5) = N.z();

        b(i) = N.x()*D.x() + N.y()*D.y() + N.z()*D.z() - N.x()*S.x() - N.y()*S.y() - N.z()*S.z();
    }

    // Solve x
    // x = (alpha beta gamma t_x t_y t_z)'T
    Eigen::MatrixXd x = ((A.transpose() * A).inverse()) * (A.transpose()) * b;

    // Compute rigid transform R and T
    Eigen::Matrix3d R;
    R.setZero();

    double sin_alpha = sin(x(0));
    double cos_alpha = cos(x(0));
    double sin_beta = sin(x(1));
    double cos_beta = cos(x(1));
    double sin_gamma = sin(x(2));
    double cos_gamma = cos(x(2));
    R(0,0) = cos_gamma * cos_beta;
    R(0,1) = -sin_gamma * cos_alpha + cos_gamma * sin_beta * sin_alpha;
    R(0,2) = sin_gamma * sin_alpha + cos_gamma * sin_beta * cos_alpha;
    R(1,0) = sin_gamma * cos_beta;
    R(1,1) = cos_gamma * cos_alpha + sin_gamma * sin_beta * sin_alpha;
    R(1,2) = -cos_gamma * sin_alpha + sin_gamma * sin_beta * cos_alpha;
    R(2,0) = -sin_beta;
    R(2,1) = cos_beta * sin_alpha;
    R(2,2) = cos_beta * cos_alpha;

    Eigen::RowVector3d T(x(3),x(4),x(5));

    return std::pair<Eigen::Matrix3d, Eigen::RowVector3d>(R,T);

}

Eigen::MatrixXd ICP::ApplyRigidTransform(Eigen::MatrixXd V_to_process, std::pair<Eigen::Matrix3d, Eigen::RowVector3d> transform){

    // #Input: V2, R, T
    // #Output: V2_Transformed

    Eigen::MatrixXd V_out;
    V_out.resize(V_to_process.rows(), V_to_process.cols());
    V_out.setZero();

    // According to the formula, p = Rq + t
    for (size_t i=0;i<V_to_process.rows();i++){
        Eigen::Vector3d row = V_to_process.row(i);
        V_out.row(i) = (transform.first * row).transpose() + transform.second;
    }

    return V_out;
}

Eigen::MatrixXd ICP::ICPOptimised(Eigen::MatrixXd V_target, Eigen::MatrixXd V_to_process, double subsample_rate){
    Eigen::MatrixXd V_subsampled = GetSubsample(V_to_process, subsample_rate);
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> correspondences = FindCorrespondences(V_target, V_subsampled);
    std::pair<double, std::pair<Eigen::Matrix3d, Eigen::RowVector3d>> transform_info = ICP::EstimateRigidTransform(correspondences.first, correspondences.second);
    return ICP::ApplyRigidTransform(V_to_process, transform_info.second);
}

//Eigen::MatrixXd ICP::ICPNormalBased(Eigen::MatrixXd V_target, Eigen::MatrixXd V_to_process){
////    Eigen::MatrixXd N = GetVertexNormal(V_target);
////    std::pair<std::pair<Eigen::MatrixXd, Eigen::MatrixXd>, Eigen::MatrixXd> correspondences = FindCorrespondencesNormalBased(V_target, V_to_process, N);
////    std::pair<Eigen::Matrix3d, Eigen::RowVector3d> transform = ICP::EstimateRigidTransformNormalBased(correspondences.first.first, correspondences.second, correspondences.first.second);
////    return ICP::ApplyRigidTransform(V_to_process, transform);
////}

double ICP::GetErrorMetric(Eigen::MatrixXd V_target, Eigen::MatrixXd V_to_process){

    double error_metric = 0.0;
    for (size_t i = 0; i < V_target.rows(); i++)
    {
        Eigen::RowVector3d V_i_target = V_target.row(i);
        Eigen::RowVector3d V_i_to_process = V_to_process.row(i);
        error_metric += (V_i_target - V_i_to_process).squaredNorm();
    }

    // Return normalised error
    return error_metric/V_target.rows();

}
