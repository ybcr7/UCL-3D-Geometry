// Global static functions

namespace MS{

    Eigen::VectorXd UniformMeanCurvature(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);
    Eigen::VectorXd UniformGaussianCurvature(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);
    Eigen::VectorXd NonUniformMeanCurvature(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);
    std::pair<Eigen::MatrixXd, Eigen::MatrixXi> Reconstruction(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);

    Eigen::SparseMatrix<double> LaplacianMatrix(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);
    Eigen::SparseMatrix<double> CotangentMatrix(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);
    Eigen::SparseMatrix<double> BarycentricArea(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);


    std::pair<Eigen::MatrixXd, Eigen::MatrixXi> ExplicitSmoothing(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);
    std::pair<Eigen::MatrixXd, Eigen::MatrixXi> ImplicitSmoothing(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);

    Eigen::MatrixXd AddNoise(Eigen::MatrixXd V_in, double noise);
}

