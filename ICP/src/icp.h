// Global static functions

namespace ICP{

    struct Transform;
    
    Eigen::MatrixXd GetSubsample(Eigen::MatrixXd V_to_process, size_t scale);
    
    Eigen::MatrixXd FindCorrespondences(Eigen::MatrixXd V_target, Eigen::MatrixXd V_to_process);
    
    Eigen::MatrixXd RejectErrors(Eigen::MatrixXd V_target, Eigen::MatrixXd V_to_process);
    
    Transform EstimateRigidTransform(Eigen::MatrixXd V_target, Eigen::MatrixXd V_source);
    
    Eigen::MatrixXd Rotate(Eigen::MatrixXd V_in, double x, double y, double z);
    
    Eigen::MatrixXd AddNoise(Eigen::MatrixXd V_in, double sd);
    
    Eigen::MatrixXd ICPBasic(Eigen::MatrixXd V_target, Eigen::MatrixXd V_to_process, size_t iteration);
    
    Eigen::MatrixXd ICPOptimised(Eigen::MatrixXd V_target, Eigen::MatrixXd V_to_process, size_t iteration);
    
    Eigen::MatrixXd ICPNormalBased(Eigen::MatrixXd V_target, Eigen::MatrixXd V_to_process, size_t iteration);
    
}

