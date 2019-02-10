// Global static functions

namespace ICP{

    std::vector<Eigen::RowVector3d> GetRigidTransform(Eigen::MatrixXd V_target, Eigen::MatrixXd V_source);
    
    Eigen::MatrixXd PointBasedICP(Eigen::MatrixXd V_target, Eigen::MatrixXd V_source, int iteration);
    
    Eigen::MatrixXd Rotate(Eigen::MatrixXd V_in, double degree);
    Eigen::MatrixXd AddNoise(Eigen::MatrixXd V_in, double sd);
    
    void PointBasedICPOptimised();
    void NormalBasedICP();
    
}
