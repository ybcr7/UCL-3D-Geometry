// Global static functions

namespace MS{

    Eigen::SparseMatrix<double> LaplacianMatrix(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);
    Eigen::SparseMatrix<double> CotangentMatrix(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);
    Eigen::SparseMatrix<double> BarycentricMassMatrix(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);
    Eigen::SparseMatrix<double> LaplacianBeltramiMatrix(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);

    Eigen::VectorXd UniformMeanCurvature(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);
    Eigen::VectorXd UniformGaussianCurvature(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);
    Eigen::VectorXd NonUniformMeanCurvature(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);

    Eigen::MatrixXd Reconstruction(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in, int k);

    Eigen::MatrixXd ExplicitSmoothing(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in, double lambda);
    Eigen::MatrixXd ImplicitSmoothing(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in, double lambda);

    Eigen::MatrixXd AddNoise(Eigen::MatrixXd V_in, double noise);
}

