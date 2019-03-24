#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <math.h>
#include <random>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/adjacency_list.h>
#include <igl/bounding_box.h>
#include "Spectra/SymEigsSolver.h"
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

    // Find connected vertex for each vertex (ordered)
    std::vector<std::vector<int>> V_adjacent;
    igl::adjacency_list(F_in, V_adjacent);

	// Find connected faces for each vertex
	std::vector<std::vector<int> > VF;
	std::vector<std::vector<int> > VFi;
	igl::vertex_triangle_adjacency(V_in.rows(), F_in, VF, VFi);

	// Compute the cotangent of founded vertex for each connected surface
    for (int i = 0; i < V_adjacent.size(); i++){

		Eigen::Vector3d V_i = V_in.row(i);
        double sum_cotangent = 0;
		
		// Find a list of faces that only connect to the current vertex
		std::vector<int> Fi_connected = VF[i];

        for (int j = 0; j < V_adjacent[i].size(); j++){

			Eigen::Vector3d V_j = V_in.row(V_adjacent[i][j]);

			std::vector<int> Vi_adjacent;

			// Find all vertices that are opposite to the current edge using face information
			for (int f = 0; f < Fi_connected.size(); f++) {
				Eigen::Vector3i F_current = F_in.row(Fi_connected[f]);

				// If any face contains this edge (i <-> j)
				if (F_current.x() == V_adjacent[i][j] || F_current.y() == V_adjacent[i][j] || F_current.z() == V_adjacent[i][j]) {

					// Construct a vertex index list for the found face
					std::vector<int> Vi_found;
					Vi_found.push_back(F_current.x());
					Vi_found.push_back(F_current.y());
					Vi_found.push_back(F_current.z());

					// Remove the edge vertex i, j and the last vertex is the one we need
					Vi_found.erase(std::remove(Vi_found.begin(), Vi_found.end(), i), Vi_found.end());
					Vi_found.erase(std::remove(Vi_found.begin(), Vi_found.end(), V_adjacent[i][j]), Vi_found.end());
					Vi_adjacent.push_back(Vi_found[0]);
				}
			}
			
			// Find the total number of founded vertices
			int Vi_adjacent_num = Vi_adjacent.size();
			double sum = 0;

			// Calculate each vertices cotangent value and accumulate it to the sum
			for (int c = 0; c < Vi_adjacent_num; c++) {
				Eigen::Vector3d V_cotagent = V_in.row(Vi_adjacent[c]);
				double cosine_gamma = (V_i - V_cotagent).dot(V_j - V_cotagent) / ((V_i - V_cotagent).norm()*(V_j - V_cotagent).norm());
				sum += 1.0 / std::tan(std::acos(cosine_gamma));
			
			}

			// Divide the sum by number of cotangent vertices
			// For a perfect manifold mesh, each edge will have exactly two connected faces (1 edge: 2 faces)
			// Then it will find two contangent values, as suggested in lecture slide, alpha and beta. Thus we need divide the sum by 2
			// Here we can also handle non-manifold (1 edge: n faces) or non-enclosed mesh (1 edge: 1 face), we need to divide the sum by n.
			sum /= Vi_adjacent_num;

			// Check cotagent matrix condition (i != j)
			if (i != V_adjacent[i][j]) {
				contangent_matrix.insert(i, V_adjacent[i][j]) = sum;
				sum_cotangent += sum;
			}
        }

        // Diagonal direction (i = j)
		contangent_matrix.insert(i, i) = -1.0 * sum_cotangent;
    }

    return contangent_matrix;

}

Eigen::SparseMatrix<double> MS::BarycentricMassMatrix(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in){

    Eigen::SparseMatrix<double> mass_matrix(V_in.rows(), V_in.rows());
	std::vector<double> area_value_list;

	// Compute area of each face using Heron's formula
	for (int i = 0; i < F_in.rows(); i++) {
		double a = (V_in.row(F_in.row(i)[0]) - V_in.row(F_in.row(i)[1])).norm();
		double b = (V_in.row(F_in.row(i)[1]) - V_in.row(F_in.row(i)[2])).norm();
		double c = (V_in.row(F_in.row(i)[2]) - V_in.row(F_in.row(i)[0])).norm();
		double s = (a + b + c) / 2;
		double A = sqrt(s*(s - a)*(s - b)*(s - c));
		area_value_list.push_back(A);
	}

    // Find connected faces for each vertex
    std::vector<std::vector<int> > VF;
    std::vector<std::vector<int> > VFi;
    igl::vertex_triangle_adjacency(V_in.rows(), F_in, VF, VFi);

    // Compute area of each connected face and sum them up
    for (int i = 0; i < VF.size(); i++){
        double total_area = 0;
        std::vector<int> faces = VF[i];
        for (int j = 0; j < faces.size(); j++){
            total_area += area_value_list[faces[j]];
        }

        // Assign the value diagonally
        mass_matrix.insert(i,i) = (total_area/3);
    }
    return mass_matrix;
}

Eigen::SparseMatrix<double> MS::LaplaceBeltramiMatrix(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in){

	Eigen::SparseMatrix<double> L_sparse(V_in.rows(), V_in.rows());

    Eigen::SparseMatrix<double> mass, cotangent;
    cotangent = CotangentMatrix(V_in, F_in);
    mass = BarycentricMassMatrix(V_in, F_in);

	// Compute the inverse of M
	Eigen::SparseMatrix<double> mass_inverse = mass.cwiseInverse();

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

Eigen::VectorXd MS::GaussianCurvature(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in){

    Eigen::VectorXd K(V_in.rows());
    Eigen::SparseMatrix<double> area = BarycentricMassMatrix(V_in, F_in);

	// Find connected faces for each vertex
	std::vector<std::vector<int> > VF;
	std::vector<std::vector<int> > VFi;
	igl::vertex_triangle_adjacency(V_in.rows(), F_in, VF, VFi);

    // Compute the Gaussian curvature based on the degree and area
    Eigen::Vector3d V_current, V_1, V_2;
    for (int i = 0; i < VF.size(); i++){
        double sum_theta = 0;
        V_current = V_in.row(i);

		// Connected vertex can be found using face information
        for (int j = 0; j < VF[i].size(); j++){

			// Construct a vertex index list for the current face
			Eigen::Vector3i F_current = F_in.row(VF[i][j]);
			std::vector<int> Vi_found;
			Vi_found.push_back(F_current.x());
			Vi_found.push_back(F_current.y());
			Vi_found.push_back(F_current.z());

			// Remove the current vertex from the current face so that we have two connected vertex
			Vi_found.erase(std::remove(Vi_found.begin(), Vi_found.end(), i), Vi_found.end());
			V_1 = V_in.row(Vi_found[0]);
			V_2 = V_in.row(Vi_found[1]);

			// Compute using the formula
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

    Eigen::SparseMatrix<double> LB_sparse = LaplaceBeltramiMatrix(V_in, F_in);
    Eigen::MatrixXd LB = LB_sparse.toDense();
    Eigen::MatrixXd cotangent_vertex = LB * V_in;

	// Compute mean curvature
    H = 0.5 * (cotangent_vertex).rowwise().norm();

    return H;
}

Eigen::MatrixXd MS::Reconstruction(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in, int k){

	Eigen::MatrixXd V_recon(V_in.rows(), V_in.cols());
	V_recon.setZero();

	Eigen::SparseMatrix<double> cotangent = CotangentMatrix(V_in, F_in);
	Eigen::SparseMatrix<double> mass = BarycentricMassMatrix(V_in, F_in);
	Eigen::SparseMatrix<double> mass_inverse_half = mass.cwiseSqrt().cwiseInverse();

	// L' = M^-0.5 * -C * M^-0.5
	Eigen::SparseMatrix<double> decomp_matrix = mass_inverse_half  * -1.0 * cotangent * mass_inverse_half;

    // Construct matrix operation object using the wrapper class SparseGenMatProd
    Spectra::SparseSymMatProd<double> operation(decomp_matrix);

    // Construct eigen solver object, requesting the largest three eigenvalues
    Spectra::SymEigsSolver< double, Spectra::SMALLEST_ALGE, Spectra::SparseSymMatProd<double> > eigen_solver(&operation, k, 10*k);

    // Initialize and compute
    eigen_solver.init();
    int nconv = eigen_solver.compute();

    // Retrieve results
    Eigen::MatrixXcd eigenvectors_complex;

    if(eigen_solver.info() == Spectra::SUCCESSFUL){
        eigenvectors_complex = eigen_solver.eigenvectors();
    }else{
        std::cout << "ERROR: No smallest eigenvectors found" << std::endl;
		return V_recon;
    }

	// Convert complex values to real values
    Eigen::MatrixXd eigenvectors(eigenvectors_complex.rows(), eigenvectors_complex.cols());

	// ei = M^-0.5 * yi
	eigenvectors = mass_inverse_half * eigenvectors_complex.real();

    for(int i=0; i<eigenvectors.cols(); i++){
		// Redefine the inner product (normalisation)
        V_recon += eigenvectors.col(i)*(V_in.transpose() * mass * eigenvectors.col(i)).transpose();
    }
	
    return V_recon;
}

Eigen::MatrixXd MS::ExplicitSmoothing(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in, double lambda, int iteration){

    Eigen::SparseMatrix<double> L = MS::LaplaceBeltramiMatrix(V_in, F_in);
    Eigen::SparseMatrix<double> I(V_in.rows(),V_in.rows());
    I.setIdentity();
    Eigen::MatrixXd V_out = V_in;

	// Compute using explicit scheme
    for (int i =0; i < iteration; i++){
        V_out = (I+lambda*L)*V_out;
    }

    return V_out;
}

Eigen::MatrixXd MS::ImplicitSmoothing(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in, double lambda, int iteration){

    Eigen::MatrixXd V_out = V_in;

    // Compute A
    Eigen::SparseMatrix<double> M = BarycentricMassMatrix(V_in, F_in);
    Eigen::SparseMatrix<double> C = CotangentMatrix(V_in, F_in);

    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> cholesky(M-lambda*C);

    // Solve the SLE
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
        double x =gaussian(rnd)*noise_scale_x/100;
        double y =gaussian(rnd)*noise_scale_y/100;
        double z =gaussian(rnd)*noise_scale_z/100;
        Eigen::RowVector3d gaussian_noise(x,y,z);
        V_out.row(i) = V_in.row(i) + gaussian_noise;
    }
    
    return V_out;
    
}