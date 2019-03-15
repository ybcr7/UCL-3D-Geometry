#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/Sparse>
#include <random>
#include <iostream>
#include <math.h>
#include <time.h>
#include <igl/boundary_loop.h>
#include <igl/bounding_box.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/adjacency_list.h>
#include <igl/fit_plane.h>
#include <igl/doublearea.h>
#include "nanoflann.hpp"
#include "ms.h"

Eigen::SparseMatrix<double> MS::UniformLaplacianMatrix(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in) {
    Eigen::SparseMatrix<double> laplacian_matrix(V_in.rows(), V_in.rows());
    std::vector<std::vector<int>> V_adjacent;
    igl::adjacency_list(F_in, V_adjacent);

    for (int i = 0; i < V_in.rows(); i++){
        int num_adjacent_vertex = V_adjacent[i].size();

        for (int j = 0; j < num_adjacent_vertex; j ++){
            laplacian_matrix.insert(i, V_adjacent[i][j]) = 1.0 / num_adjacent_vertex;
        }

        laplacian_matrix.insert(i,i) = -1.0;
    }

    return laplacian_matrix;
}

Eigen::SparseMatrix<double> MS::NonUniformLaplacianMatrix(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in) {
    Eigen::SparseMatrix<double> laplacian_matrix(V_in.rows(), V_in.rows());


//    for (int i = 0; i < V_in.rows(); i++){
//        int num_adjacent_vertex = V_adjacent.size();
//
//        for (int j = 0; j < num_adjacent_vertex; j ++){
//            laplacian_matrix.insert(i, V_adjacent[i][j]) = 1.0 / num_adjacent_vertex;
//        }
//
//        laplacian_matrix.insert(i,i) = -1.0;
//
//    }
//
//    return laplacian_matrix;

}


Eigen::SparseMatrix<double> MS::BarycentricArea(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in){
    Eigen::SparseMatrix<double> area_sparse;

    Eigen::MatrixXd area_each_vertex(V_in.rows(),1);

    // Compute the area for each face in the mesh
    Eigen::MatrixXd area_value_list;
    igl::doublearea(V_in, F_in, area_value_list);

    // Find connected faces for each vertex
    std::vector<std::vector<int> > VF;
    std::vector<std::vector<int> > VFi;
    igl::vertex_triangle_adjacency(V_in.rows(), F_in, VF, VFi);
    for (int i = 0; i < VF.size(); i++){
        double total_area = 0;
        std::vector<int> faces = VF[i];
        for (int j = 0; j < faces.size(); j++){
            total_area += area_value_list.row(faces[j])[0];
        }
        area_each_vertex.row(i)[0] = total_area/3;
    }

    // Convert the column to diagonal
    Eigen::MatrixXd area_dense(area_each_vertex.asDiagonal());

    // Convert the dense matrix to sparse
    area_sparse = area_dense.sparseView();

    return area_sparse;
}



Eigen::VectorXd MS::UniformMeanCurvature(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in){

    Eigen::VectorXd H(V_in.rows());
    H.setZero();

    Eigen::MatrixXd laplacian_dense = UniformLaplacianMatrix(V_in,F_in).toDense();
    Eigen::MatrixXd laplacian_vertex = laplacian_dense * V_in;

    H = 0.5 * (laplacian_vertex).rowwise().norm();
    return H;
}

Eigen::VectorXd MS::UniformGaussianCurvature(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in){

}

Eigen::VectorXd MS::NonUniformMeanCurvature(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in){

    Eigen::VectorXd H(V_in.rows());
    H.setZero();

    Eigen::SparseMatrix<double> area_inverse,laplacian;

    laplacian =  UniformLaplacianMatrix(V_in, F_in);
    area_inverse = BarycentricArea(V_in, F_in);

    Eigen::SparseMatrix<double> AreaInv(V_in.rows(), V_in.rows());

    for (int i = 0; i < V_in.rows(); i++)
    {
        AreaInv.insert(i, i) = 1 / area_inverse.coeff(i, i);
    }

    Eigen::SparseMatrix<double> sparse_matrix(V_in.rows(), V_in.rows());
    sparse_matrix = AreaInv * laplacian;
    //cout << "sparse_matrix " << sparse_matrix.rows() << endl;
    Eigen::MatrixXd dsM = sparse_matrix.toDense();
    Eigen::MatrixXd laplacian_vertex = dsM * V_in;

    H = 0.5 * (laplacian_vertex).rowwise().norm();
    return H;

}


Eigen::MatrixXd MS::AddNoise(Eigen::MatrixXd V_in, double noise){
    
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
    std::normal_distribution<double> gaussian(0.0, noise);
    
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
