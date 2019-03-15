// Custom Data Structure
#include <eigen/Eigen/src/Core/Matrix.h>

namespace CDS{
    struct RenderingData{
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        Eigen::MatrixXd C;
    };
}