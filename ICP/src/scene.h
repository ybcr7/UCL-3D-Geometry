// Scene manager for loading, allocating tasks and displaying meshes
class Scene
{
public:
    Scene(igl::opengl::glfw::Viewer& refViewer);
    ~Scene();
    
    // Task 1
    void Point2PointAlign();
    
    // Task 2 & 3
    void RotateMeshWithNoise(double x, double y, double z, double sd);
    
    // Task 4
    void Point2PointAlignOptimised();
    
    // Task 5
    void MuiltMeshAlign();
    
    // Task 6
    void Point2PlaneAlign();
    
    // Utility
    void Initialise();
    void DisplayMeshes();
    
    void SetIteration(int i);
    
private:
    
    igl::opengl::glfw::Viewer& viewer;
    
    Eigen::MatrixXd V1, V2, V3, V4, V5;
    Eigen::MatrixXi F1, F2, F3, F4, F5;
    
    std::vector<Eigen::MatrixXd> Vs;
    std::vector<Eigen::MatrixXi> Fs;
    
    int iteration;
    
};
