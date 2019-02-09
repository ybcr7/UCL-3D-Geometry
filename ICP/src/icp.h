// Global static functions

namespace ICP{

    Eigen::MatrixXd Rotate(Eigen::MatrixXd V_in, double degree);
    Eigen::MatrixXd AddNoise(Eigen::MatrixXd V_in, double sd);
    void PointBasedICP();
    void PointBasedICPOptimised();
    void NormalBasedICP();
    
}
