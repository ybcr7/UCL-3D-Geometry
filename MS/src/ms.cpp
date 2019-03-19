#define _USE_MATH_DEFINES
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <random>
#include <iostream>
#include <math.h>
#include <time.h>
#include <igl/bounding_box.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/adjacency_list.h>
#include <igl/doublearea.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/cotmatrix_entries.h>
#include <igl/is_edge_manifold.h>
#include <igl/is_symmetric.h>
#include <igl/eigs.h>
#include "Spectra/GenEigsSolver.h"
#include "Spectra/SymEigsSolver.h"
#include "Spectra/MatOp/SparseGenMatProd.h"
#include "Spectra/MatOp/SparseSymMatProd.h"
#include "ms.h"

Eigen::SparseMatrix<double> MS::LaplacianMatrix(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in) {

    Eigen::SparseMatrix<double> laplacian_matrix(V_in.rows(), V_in.rows());

    // Find connected vertex for each vertex
    std::vector<std::vector<int>> V_adjacent;
    igl::adjacency_list(F_in, V_adjacent);

    // Construct the Laplacian matrix based on the number of neighbours
    for (int i = 0; i < V_in.rows(); i++){
        int num_adjacent_vertex = V_adjacent[i].size();

        for (int j = 0; j < num_adjacent_vertex; j ++){
            laplacian_matrix.insert(i, V_adjacent[i][j]) = 1.0 / num_adjacent_vertex;
        }

        // Assign -1 to diagonal direction
        laplacian_matrix.insert(i,i) = -1.0;
    }

    return laplacian_matrix;
}

Eigen::SparseMatrix<double> MS::CotangentMatrix(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in){

    Eigen::SparseMatrix<double> contangent_matrix(V_in.rows(), V_in.rows());

    if (igl::is_edge_manifold(F_in)){
		// Manifold
		std::cout << "Manifold Mesh - Using MS::CotangentMatrix()" << std::endl;
        
        // Find connected vertex for each vertex (ordered)
        std::vector<std::vector<int>> V_adjacent;
        igl::adjacency_list(F_in, V_adjacent, true);

		// Find connected faces for each vertex
		std::vector<std::vector<int> > VF;
		std::vector<std::vector<int> > VFi;
		igl::vertex_triangle_adjacency(V_in.rows(), F_in, VF, VFi);

        // Compute the cotangent of alpha and beta and construct the Laplacian matrix
        Eigen::Vector3d V_current, V_1, V_2, V_3;
        for (int i = 0; i < V_adjacent.size(); i++){

            double sum_cotangent = 0;
            V_current = V_in.row(i);

			// Find a list of faces that only connect to the current vertex
			std::vector<int> Fi_connected = VF[i];

            // Find three continuous vertex in order to compute the degrees
            for (int j = 0; j < V_adjacent[i].size(); j++){

				std::vector<int> Vi_adjacent_pair;

				// Find two connected faces with and target vertex
				for (int f = 0; f < Fi_connected.size(); f++) {
					Eigen::Vector3i ff = F_in.row(Fi_connected[f]);
					if (ff.x() == V_adjacent[i][j] || ff.y() == V_adjacent[i][j] || ff.z() == V_adjacent[i][j]) {
						//std::cout << "Found" << std::endl;
						std::vector<int> V_found;
						V_found.push_back(ff.x());
						V_found.push_back(ff.y());
						V_found.push_back(ff.z());
						V_found.erase(std::remove(V_found.begin(), V_found.end(), i), V_found.end());
						V_found.erase(std::remove(V_found.begin(), V_found.end(), V_adjacent[i][j]), V_found.end());
						Vi_adjacent_pair.push_back(V_found[0]);
					}
				}

				if (Vi_adjacent_pair.size() != 2) {

					std::cout << "ERROR" << std::endl;
				}

				V_1 = V_in.row(Vi_adjacent_pair[0]);
				V_2 = V_in.row(V_adjacent[i][j]);
				V_3 = V_in.row(Vi_adjacent_pair[1]);

                double cosine_alpha = (V_current - V_1).dot(V_2 - V_1)/((V_current - V_1).norm()*(V_2 - V_1).norm());
                double cosine_beta = (V_current - V_3).dot(V_2 - V_3)/((V_current - V_3).norm()*(V_2 - V_3).norm());
                double sum = 0.5 * (1.0/std::tan(std::acos(cosine_alpha)) + 1.0/std::tan(std::acos(cosine_beta)));

                // If I != J
                if (V_in.row(i) != V_in.row(V_adjacent[i][j])){
                    contangent_matrix.insert(i, V_adjacent[i][j]) = sum;
                    sum_cotangent += sum;
                }

            }
            // Diagonal direction (I = J)
            contangent_matrix.insert(i,i) = -1.0 * sum_cotangent;
        }
		std::cout << "What" << std::endl;

        return contangent_matrix;

    }else{
        // Non-Manifold
		std::cout << "Non-Manifold Mesh - Using igl::cotmatrix()" << std::endl;
		igl::cotmatrix(V_in, F_in, contangent_matrix);
        return contangent_matrix;
    };

}

Eigen::SparseMatrix<double> MS::BarycentricMassMatrix(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in){

    Eigen::SparseMatrix<double> mass_matrix(V_in.rows(), V_in.rows());

    // Compute the area for each face in the mesh
    Eigen::MatrixXd area_value_list;

    // This will compute twice the area for each face thus we need to half it
    igl::doublearea(V_in, F_in, area_value_list);
    area_value_list *= 0.5;

    // Find connected faces for each vertex
    std::vector<std::vector<int> > VF;
    std::vector<std::vector<int> > VFi;
    igl::vertex_triangle_adjacency(V_in.rows(), F_in, VF, VFi);

	std::cout << VF.size() << std::endl;
	std::cout << VFi.size() << std::endl;

    // Compute area of each connected face and sum them up
    for (int i = 0; i < VF.size(); i++){
        double total_area = 0;
        std::vector<int> faces = VF[i];
        for (int j = 0; j < faces.size(); j++){
            total_area += area_value_list.row(faces[j])[0];
        }

        // Assign the value diagonally
        mass_matrix.insert(i,i) = (total_area/3);
    }
    return mass_matrix;
}

Eigen::SparseMatrix<double> MS::LaplacianBeltramiMatrix(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in){

    Eigen::SparseMatrix<double> mass, cotangent;
    cotangent = CotangentMatrix(V_in, F_in);
    mass = BarycentricMassMatrix(V_in, F_in);

    // For testing purpose: To show that we can obtain the exactly same result obtained from our implementations versus built-in functions
//    Eigen::SparseMatrix<double> mass_built_in, cotangent_built_in;
//    igl::massmatrix(V_in, F_in, igl::MASSMATRIX_TYPE_BARYCENTRIC, mass_built_in);
//    igl::cotmatrix(V_in, F_in, cotangent_built_in);
//    std::cout << "Result(M): Custom vs Built-in Function" << std::endl;
//    std::cout << mass.row(0) << std::endl;
//    std::cout << "-----------------------------------------" << std::endl;
//    std::cout << mass_built_in.row(0) << std::endl;
//    std::cout << "Result(C): Custom vs Built-in Function" << std::endl;
//    Eigen::VectorXd cot_row1 = cotangent.row(0);
//    Eigen::VectorXd cot_row2 = cotangent_built_in.row(0);
//    for (int i =0; i < cot_row1.size(); i ++){
//        if (cot_row1[i]!=0){
//            std::cout << cot_row1[i] << std::endl;
//        }
//    }
//    std::cout << "-----------------------------------------" << std::endl;
//    for (int i =0; i < cot_row2.size(); i ++){
//        if (cot_row2[i]!=0){
//            std::cout << cot_row2[i] << std::endl;
//        }
//    }

	Eigen::SparseMatrix<double> mass_inverse = mass.cwiseInverse();

    Eigen::SparseMatrix<double> L_sparse(V_in.rows(), V_in.rows());
    L_sparse = mass_inverse * cotangent;

    return L_sparse;
}


Eigen::VectorXd MS::UniformMeanCurvature(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in){

    Eigen::VectorXd H(V_in.rows());
    H.setZero();

    Eigen::MatrixXd laplacian_dense = LaplacianMatrix(V_in,F_in).toDense();
    Eigen::MatrixXd laplacian_vertex = laplacian_dense * V_in;

	// Compute mean curvature
    H = 0.5 * (laplacian_vertex).rowwise().norm();

    return H;
}

Eigen::VectorXd MS::UniformGaussianCurvature(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in){

    Eigen::VectorXd K(V_in.rows());
    Eigen::SparseMatrix<double> area = BarycentricMassMatrix(V_in, F_in);


	// Find connected faces for each vertex
	std::vector<std::vector<int> > VF;
	std::vector<std::vector<int> > VFi;
	igl::vertex_triangle_adjacency(V_in.rows(), F_in, VF, VFi);

    // Find connected vertex for each vertex (ordered)
    std::vector<std::vector<int>> V_adjacent;
    igl::adjacency_list(F_in, V_adjacent, true);

    // Compute the Gaussian curvature based on the degree and area
    Eigen::Vector3d V_current, V_1, V_2;
    for (int i = 0; i < V_adjacent.size(); i++){
        double sum_theta = 0;
        V_current = V_in.row(i);
        for (int j = 0; j < V_adjacent[i].size(); j++){

            // Find two continuous vertex in order to compute the degree between two vectors
            if (j + 1 == V_adjacent[i].size()){
                V_1 = V_in.row(V_adjacent[i][j]);
                V_2 = V_in.row(V_adjacent[i][0]);
            }else{
                V_1 = V_in.row(V_adjacent[i][j]);
                V_2 = V_in.row(V_adjacent[i][j+1]);
            }

            double theta = std::acos((V_1 - V_current).dot(V_2 - V_current)/((V_1 - V_current).norm() * (V_2 - V_current).norm()));
            sum_theta += theta;
        }

        // Map the Gaussian curvature to each vertex
        K[i] = (2*M_PI-sum_theta)/area.coeff(i,i);
    }

    return K;
}

Eigen::VectorXd MS::NonUniformMeanCurvature(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in){

    Eigen::VectorXd H(V_in.rows());
    H.setZero();

    Eigen::SparseMatrix<double> LB_sparse = LaplacianBeltramiMatrix(V_in, F_in);

    Eigen::MatrixXd LB = LB_sparse.toDense();
    Eigen::MatrixXd cotangent_vertex = LB * V_in;

	// Compute mean curvature
    H = 0.5 * (cotangent_vertex).rowwise().norm();

    return H;
}

Eigen::MatrixXd MS::Reconstruction(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in, int k){

	//Eigen::SparseMatrix<double> cotangent = CotangentMatrix(V_in, F_in);

	Eigen::SparseMatrix<double> cotangent;
	igl::cotmatrix(V_in, F_in, cotangent);

	//Eigen::SparseMatrix<double> mass2;
	//igl::massmatrix(V_in, F_in, igl::MASSMATRIX_TYPE_BARYCENTRIC, mass2);
	Eigen::SparseMatrix<double> mass = BarycentricMassMatrix(V_in, F_in);
	
	Eigen::SparseMatrix<double> mass_inverse = mass.cwiseInverse();
	Eigen::SparseMatrix<double> mass_inverse_half = mass.cwiseSqrt().cwiseInverse();

	Eigen::SparseMatrix<double> decomp_matrix = mass_inverse_half  * -1.0 * cotangent * mass_inverse_half;

    // Construct matrix operation object using the wrapper class SparseGenMatProd
    Spectra::SparseSymMatProd<double> operation(decomp_matrix);

    // Construct eigen solver object, requesting the largest three eigenvalues
    Spectra::SymEigsSolver< double, Spectra::SMALLEST_ALGE, Spectra::SparseSymMatProd<double> > eigen_solver(&operation, k, 10*k);

    // Initialize and compute
    eigen_solver.init();
    int nconv = eigen_solver.compute();
    // Retrieve results
    Eigen::VectorXcd eigenvalues_complex;
    Eigen::MatrixXcd eigenvectors_complex;

    if(eigen_solver.info() == Spectra::SUCCESSFUL){
        eigenvalues_complex = eigen_solver.eigenvalues();
        eigenvectors_complex = eigen_solver.eigenvectors();
    }else{
        std::cout << "ERROR: No smallest eigenvectors found" << std::endl;
		return V_in;
    }

	// Convert complex values to real values
    Eigen::MatrixXd eigenvectors(eigenvectors_complex.rows(), eigenvectors_complex.cols());

	eigenvectors = eigenvectors_complex.real();
	eigenvectors = mass_inverse_half * eigenvectors_complex.real();

    Eigen::MatrixXd V_recon(V_in.rows(),3);
    V_recon.setZero();

    for(int i=0; i<eigenvectors.cols(); i++){
        V_recon += eigenvectors.col(i)*(V_in.transpose() * mass * eigenvectors.col(i)).transpose();

		// Equivalent to:
		//V_recon.col(0) += (V_in.col(0).transpose() * eigenvectors.col(i))* eigenvectors.col(i);
		//V_recon.col(1) += (V_in.col(1).transpose() * eigenvectors.col(i))* eigenvectors.col(i);
		//V_recon.col(2) += (V_in.col(2).transpose() * eigenvectors.col(i))* eigenvectors.col(i);
    }
	
    return V_recon;
}

Eigen::MatrixXd MS::ExplicitSmoothing(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in, double lambda, int iteration){

    //Eigen::SparseMatrix<double> L = MS::LaplacianMatrix(V_in, F_in);
    Eigen::SparseMatrix<double> L = MS::LaplacianBeltramiMatrix(V_in, F_in);
    Eigen::SparseMatrix<double> I(V_in.rows(),V_in.rows());
    I.setIdentity();
    Eigen::MatrixXd V_out = V_in;

    for (int i =0; i < iteration; i++){
        V_out = (I+lambda*L)*V_out;
    }

    return V_out;
}

Eigen::MatrixXd MS::ImplicitSmoothing(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in, double lambda, int iteration){

    Eigen::MatrixXd V_out = V_in;

    // Compute A
    Eigen::SparseMatrix<double> L = LaplacianBeltramiMatrix(V_in, F_in);
    Eigen::SparseMatrix<double> M = BarycentricMassMatrix(V_in, F_in);
    Eigen::SparseMatrix<double> C = CotangentMatrix(V_in, F_in);

    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> cholesky(M-lambda*C);
    //Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> cholesky(M-lambda*M*L);

    // Compute b
    for (int i =0; i< iteration; i++){
        V_out = cholesky.solve(M*V_out);
    }

    return V_out;
}

Eigen::MatrixXd MS::AddNoise(Eigen::MatrixXd V_in, double sd){
    
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
        Eigen::RowVector3d gaussian_noise(x,y,z);
        V_out.row(i) = V_in.row(i) + gaussian_noise;
    }
    
    return V_out;
    
}

void MS::test(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in){

    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> C_out;

//    igl::cotmatrix_entries(V_in,F_in,C_out);

//    std::cout << F_in.rows() << std::endl;
//
//    std::cout << C_out.rows() << std::endl;
//
//    std::cout << C_out.row(1) << std::endl;
    //std::cout << C_out << std::endl;

}



