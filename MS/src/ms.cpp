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
#include <igl/cotmatrix.h>
#include "nanoflann.hpp"
#include "ms.h"

Eigen::SparseMatrix<double> MS::LaplacianMatrix(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in) {
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

Eigen::SparseMatrix<double> MS::CotangentMatrix(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in){
    Eigen::SparseMatrix<double> contangent_sparse;

    Eigen::Vector3d cotangent_sum_each_vertex(V_in.rows());
    Eigen::SparseMatrix<double> area = BarycentricArea(V_in, F_in);

    std::vector<std::vector<int>> V_adjacent;
    igl::adjacency_list(F_in, V_adjacent, true);

    Eigen::Vector3d V_current, V_1, V_2, V_3;

    for (int i = 0; i < V_adjacent.size(); i++){
        double sum_cotangent = 0;
        V_current = V_in.row(i);
        for (int j = 0; j < V_adjacent[i].size(); j++){
            if (j == 0){
                V_1 = V_in.row(V_adjacent[i][V_adjacent[i].back()]);
                V_2 = V_in.row(V_adjacent[i][j]);
                V_3 = V_in.row(V_adjacent[i][j+1]);
            }
            else if (j == V_adjacent[i].size()){
                V_1 = V_in.row(V_adjacent[i][j-1]);
                V_1 = V_in.row(V_adjacent[i][j]);
                V_2 = V_in.row(V_adjacent[i][0]);
            }else{
                V_1 = V_in.row(V_adjacent[i][j-1]);
                V_2 = V_in.row(V_adjacent[i][j]);
                V_3 = V_in.row(V_adjacent[i][j+1]);
            }

            double cosine_alpha = (V_current - V_1).dot(V_2 - V_1)/((V_current - V_1).norm()*(V_2 - V_1).norm());
            double cosine_beta = (V_current - V_3).dot(V_2 - V_3)/((V_current - V_3).norm()*(V_2 - V_3).norm());
            double sum = cosine_alpha/std::sin(std::acos(cosine_alpha)) + cosine_beta/std::sin(std::acos(cosine_beta));

            sum_cotangent += sum;
        }
        cotangent_sum_each_vertex[i] = sum_cotangent;
    }

    // Convert the column to diagonal
    Eigen::MatrixXd cotangent_dense(cotangent_sum_each_vertex.asDiagonal());

    // Convert the dense matrix to sparse
    contangent_sparse = cotangent_dense.sparseView();

    return contangent_sparse;

}

Eigen::SparseMatrix<double> MS::BarycentricArea(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in){
    Eigen::SparseMatrix<double> area_sparse;

    Eigen::VectorXd area_each_vertex(V_in.rows());

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
        area_each_vertex[i] = total_area/3;
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

    Eigen::MatrixXd laplacian_dense = LaplacianMatrix(V_in,F_in).toDense();
    Eigen::MatrixXd laplacian_vertex = laplacian_dense * V_in;

    H = 0.5 * (laplacian_vertex).rowwise().norm();
    return H;
}

Eigen::VectorXd MS::UniformGaussianCurvature(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in){

    Eigen::VectorXd k_each_vertex(V_in.rows());
    Eigen::SparseMatrix<double> area = BarycentricArea(V_in, F_in);

    std::vector<std::vector<int>> V_adjacent;
    igl::adjacency_list(F_in, V_adjacent, true);

    Eigen::Vector3d V_current, V_1, V_2;

    for (int i = 0; i < V_adjacent.size(); i++){
        double sum_theta = 0;
        V_current = V_in.row(i);
        for (int j = 0; j < V_adjacent[i].size(); j++){

            if (j == V_adjacent[i].size()){
                V_1 = V_in.row(V_adjacent[i][j]);
                V_2 = V_in.row(V_adjacent[i][V_adjacent[i].front()]);
            }else{
                V_1 = V_in.row(V_adjacent[i][j]);
                V_2 = V_in.row(V_adjacent[i][j+1]);
            }

            double theta = std::acos((V_1 - V_current).dot(V_2 - V_current)/((V_1 - V_current).norm() * (V_2 - V_current).norm()));
            sum_theta += theta;
        }
        k_each_vertex[i] = (2*M_PI-sum_theta)/area.coeff(i,i);
    }

    return k_each_vertex;
}

Eigen::VectorXd MS::NonUniformMeanCurvature(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in){

    Eigen::VectorXd H(V_in.rows());
    H.setZero();

    Eigen::SparseMatrix<double> area,laplacian;

    laplacian = CotangentMatrix(V_in, F_in);
    //igl::cotmatrix(V_in, F_in, laplacian);
    //laplacian =  LaplacianMatrix(V_in, F_in);
    area = BarycentricArea(V_in, F_in);

    Eigen::SparseMatrix<double> area_inverse(V_in.rows(), V_in.rows());

    for (int i = 0; i < V_in.rows(); i++)
    {
        area_inverse.insert(i, i) = 1 / area.coeff(i, i);
    }

    Eigen::SparseMatrix<double> sparse_matrix(V_in.rows(), V_in.rows());
    sparse_matrix = area_inverse * laplacian;
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

