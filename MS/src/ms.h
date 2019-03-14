// Global static functions

namespace MS{

    std::pair<Eigen::MatrixXd, Eigen::MatrixXi> UniformMeanCurvature(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);
    std::pair<Eigen::MatrixXd, Eigen::MatrixXi> UniformGaussianCurvature(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);
    std::pair<Eigen::MatrixXd, Eigen::MatrixXi> NonUniformMeanCurvature(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);
    std::pair<Eigen::MatrixXd, Eigen::MatrixXi> Reconstruction(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);

    std::pair<Eigen::MatrixXd, Eigen::MatrixXi> ExplicitSmoothing(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);
    std::pair<Eigen::MatrixXd, Eigen::MatrixXi> ImplicitSmoothing(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);

    Eigen::MatrixXd AddNoise(Eigen::MatrixXd V_in, double sd);
}

