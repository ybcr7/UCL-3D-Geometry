#include <Eigen/Dense>
#include <random>
#include <iostream>
#include "nanoflann.hpp"
#include "icp.h"

std::vector<Eigen::RowVector3d> GetRigidTransform(Eigen::MatrixXd V_target, Eigen::MatrixXd V_source){
    
    
    
    
    
    
    
    
}


Eigen::MatrixXd ICP::PointBasedICP(Eigen::MatrixXd V_target, Eigen::MatrixXd V_to_process, int iteration){
    
    Eigen::MatrixXd V_out;
    V_out.resize(V_to_process.rows(), V_to_process.cols());
    V_out.setZero();

    const size_t num_result = 1;
    const size_t max_leaf = 10;
    
    nanoflann::KDTreeEigenMatrixAdaptor<Eigen::MatrixXd> kd_tree_index(V_target, max_leaf);
    kd_tree_index.index->buildIndex();
    
    for (size_t v=0; v<V_out.rows(); v++){
        
    // Pick the current vertex for query
    Eigen::RowVector3d query_vertex = V_out.row(v);
    
    // Create a query object
    std::vector<size_t> indexes(num_result);
    std::vector<double> dists_sqr(num_result);
    
    // Find the closest 1 vertex
    nanoflann::KNNResultSet<double> result(num_result);
    result.init(indexes.data(), dists_sqr.data());
    kd_tree_index.index->findNeighbors(result, query_vertex.data(), nanoflann::SearchParams(max_leaf));
    
    V_out.row(v) = V_to_process.row(indexes[0]);

    }
    
    std::cout << "Aligned" << std::endl;
    
    return V_out;
}


Eigen::MatrixXd ICP::Rotate(Eigen::MatrixXd V_in, double degree){
    
    // Initialise
    Eigen::MatrixXd V_out;
    V_out.resize(V_in.rows(), V_in.cols());
    V_out.setZero();
    
    // Construct rotation matrix
    Eigen::Matrix3d R;
    R << Eigen::AngleAxisd(degree * M_PI/180, Eigen::Vector3d(0,0,1)).toRotationMatrix();
    
    // Rotate the vertex
    V_out = V_in * R;
    
    return V_out;
    
}

Eigen::MatrixXd ICP::AddNoise(Eigen::MatrixXd V_in, double sd){
    
    // Initialise
    Eigen::MatrixXd V_out;
    V_out.resize(V_in.rows(), V_in.cols());
    V_out.setZero();
    
    // Random number generation
    std::default_random_engine rnd;
    std::normal_distribution<double> gaussian(0.0, sd);
    
    // Add noise to the vertex
    for (int i=0; i<V_out.rows(); i++){
        double x =gaussian(rnd)/100;
        double y =gaussian(rnd)/100;
        double z =gaussian(rnd)/100;
        Eigen::RowVector3d noise(x,y,z);
        V_out.row(i) = V_in.row(i) + noise;
    }
    
    return V_out;
    
}
