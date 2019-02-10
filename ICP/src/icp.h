// Global static functions

namespace ICP{

    std::vector<Eigen::RowVector3d> EstimateRigidTransform(Eigen::MatrixXd V_target, Eigen::MatrixXd V_source);
    
    Eigen::MatrixXd FindCorrespondences(Eigen::MatrixXd V_target, Eigen::MatrixXd V_to_process);
    
    Eigen::MatrixXd Rotate(Eigen::MatrixXd V_in, double degree);
    Eigen::MatrixXd AddNoise(Eigen::MatrixXd V_in, double sd);
    
    void PointBasedICPOptimised();
    void NormalBasedICP();
    
    Eigen::MatrixXd ICPBasic(Eigen::MatrixXd V_target, Eigen::MatrixXd V_to_process, size_t iteration);
    
    Eigen::MatrixXd ICPOptimised(Eigen::MatrixXd V_target, Eigen::MatrixXd V_to_process, size_t iteration);
    
    Eigen::MatrixXd ICPNormalBased(Eigen::MatrixXd V_target, Eigen::MatrixXd V_to_process, size_t iteration);
}
