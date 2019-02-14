#include <Eigen/Dense>
#include <Eigen/SVD>
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
    
    Transform transform;
    
    // Define Barycenters
    // p_bar = 1/m sigma p_i => average
    Eigen::RowVector3d p_bar = V_to_process.colwise().mean();
    Eigen::RowVector3d q_bar = V_matched.colwise().mean();
    
    Eigen::MatrixXd p_hat = V_to_process.rowwise() - p_bar;
    Eigen::MatrixXd q_hat = V_matched.rowwise() - q_bar;
    
    // Construct A
    Eigen::Matrix3d A;
    A.setZero();
    
    for (size_t i=0; i<V_to_process.rows(); i++){
        
        Eigen::Vector3d q_i = q_hat.row(i);
        Eigen::Vector3d p_i = p_hat.row(i);
        
        A += q_i * p_i.transpose();
    }
    
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU| Eigen::ComputeThinV);
    Eigen::MatrixXd R = svd.matrixV() * svd.matrixU().transpose();
    Eigen::RowVector3d T = p_bar - R * q_bar;

    transform.R = R;
    transform.T = T;
    
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

Eigen::MatrixXd ICP::ICPBasic(Eigen::MatrixXd V_target, Eigen::MatrixXd V_to_process){
    
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


//Eigen::MatrixXd ICP::ICPOptimised(Eigen::MatrixXd V_target, Eigen::MatrixXd V_to_process, size_t scale){
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
