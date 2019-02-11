#include <Eigen/Dense>
#include <random>
#include <iostream>
#include "nanoflann.hpp"
#include "icp.h"

struct ICP::Transform{
    Eigen::Matrix3d R;
    Eigen::RowVector3d T;
    Eigen::RowVector3d S;
};

ICP::Transform ICP::EstimateRigidTransform(Eigen::MatrixXd V_to_process, Eigen::MatrixXd V_matched){
    
    // The rigid transform can be estimiated from min(R,t) Sigma_i ||p_i - R*q_i - t||^2
    // t = p_bar - R*q_bar
    // R can be estimated from min(R) Sigma_i ||p_hat_i - R* q_hat_i||^2
    
    //std::cout << V_target.row(0) << std::endl;
    //std::cout << V_matched.row(0) << std::endl;
    
    Transform transform;
    
    // Define Barycenters
    
//    std::cout << V_to_process.rows() << std::endl;
//    std::cout << V_matched.rows() << std::endl;
    
    // p_bar = 1/m sigma p_i => average
    Eigen::RowVector3d p_bar = V_to_process.colwise().mean();
    Eigen::RowVector3d q_bar = V_matched.colwise().mean();
    
    //std::cout << p_bar << std::endl;
    //std::cout << q_bar << std::endl;
    
    Eigen::MatrixXd p_hat = V_to_process.rowwise() - p_bar;
    Eigen::MatrixXd q_hat = V_matched.rowwise() - q_bar;
    
    //std::cout << p_hat.row(0) << std::endl;
    //std::cout << q_hat.row(0) << std::endl;
    
    //std::cout << q_hat.row(0)*(p_hat.row(0).transpose()) << std::endl;
    
    // Construct A
    
    Eigen::Matrix3d A;
    A.setZero();
    
    for (size_t i=0; i<V_to_process.rows(); i++){
        
        Eigen::Vector3d q_i = q_hat.row(i);
        Eigen::Vector3d p_i = p_hat.row(i);
        
        A += q_i * p_i.transpose();
    }
    
    //std::cout << A << std::endl;
    
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::DecompositionOptions::ComputeThinU| Eigen::DecompositionOptions::ComputeThinV);
    Eigen::MatrixXd R = svd.matrixV() * svd.matrixU().transpose();
    Eigen::RowVector3d T = p_bar - R * q_bar;

    transform.R = R;
    transform.T = T;

    std::cout << transform.R << std::endl;
    std::cout << transform.T << std::endl;
    
    return transform;

}

Eigen::MatrixXd RejectErrors(Eigen::MatrixXd V_target, Eigen::MatrixXd V_to_process){


}

Eigen::MatrixXd ICP::FindCorrespondences(Eigen::MatrixXd V_target, Eigen::MatrixXd V_to_process){
    
    // Initialise output matrix
    Eigen::MatrixXd V_out;
    V_out.resize(V_to_process.rows(), V_to_process.cols());
    V_out.setZero();

    const size_t num_result = 1;
    const size_t max_leaf = 10;
    
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

    }
    
    return V_out;
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
    
    // Random number generation
    std::default_random_engine rnd;
    std::normal_distribution<double> gaussian(0.0, sd);
    
    // Add noise to the vertex
    for (int i=0; i<V_out.rows(); i++){
        double x =gaussian(rnd)/1000;
        double y =gaussian(rnd)/1000;
        double z =gaussian(rnd)/1000;
        Eigen::RowVector3d noise(x,y,z);
        V_out.row(i) = V_in.row(i) + noise;
    }
    
    return V_out;
    
}

Eigen::MatrixXd ICP::ICPBasic(Eigen::MatrixXd V_target, Eigen::MatrixXd V_to_process, size_t iteration){
    
    Eigen::MatrixXd V_out;
    V_out.resize(V_to_process.rows(), V_to_process.cols());
    V_out.setZero();
    
    Eigen::MatrixXd V_matched = FindCorrespondences(V_target, V_to_process);
    Transform transform = EstimateRigidTransform(V_to_process, V_matched);
    
    for (size_t i=0;i<V_to_process.rows();i++){
        
        Eigen::Vector3d row = V_to_process.row(i);
        V_out.row(i) = (row.transpose() * transform.R - transform.T).transpose();
        
    }
    
    return V_out;
    
}


//Eigen::MatrixXd ICP::ICPOptimised(Eigen::MatrixXd V_target, Eigen::MatrixXd V_to_process, size_t scale, size_t iteration){
//
//    Eigen::MatrixXd V_processed;
//
//    V_processed = GetSubsample(V_to_process, scale);
//
//    V_processed = FindCorrespondences(V_target, V_processed);
//
//    EstimateRigidTransform(V_target, V_processed);
//
//    return V_processed;
//
//}
