namespace MS{

    Eigen::SparseMatrix<double> LaplacianMatrix(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);
    Eigen::SparseMatrix<double> CotangentMatrix(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);
    Eigen::SparseMatrix<double> BarycentricMassMatrix(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);
    Eigen::SparseMatrix<double> LaplaceBeltramiMatrix(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);

    Eigen::VectorXd UniformMeanCurvature(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);
    Eigen::VectorXd GaussianCurvature(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);
    Eigen::VectorXd NonUniformMeanCurvature(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);
    Eigen::MatrixXd Reconstruction(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in, int k);

    Eigen::MatrixXd ExplicitSmoothing(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in, double lambda, int iteration);
    Eigen::MatrixXd ImplicitSmoothing(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in, double lambda, int iteration);
    Eigen::MatrixXd AddNoise(Eigen::MatrixXd V_in, double noise);

	Eigen::SparseMatrix<double> CotangentMatrix2(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);
	Eigen::SparseMatrix<double> CotangentMatrix3(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);
}

