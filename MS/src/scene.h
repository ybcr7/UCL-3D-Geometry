// Scene manager for loading, allocating tasks and displaying data
class Scene
{
public:
    Scene(igl::opengl::glfw::Viewer& refViewer);
    ~Scene();
    
    // Part A: Discrete Curvature and Spectral Meshes
    void Discretisation(int mode);
    void Reconstruction();

    // Part B: Laplacian Mesh Smoothing
    void Smoothing(int mode);
    void AddNoise();

    // Utility
    void Initialise(std::string filename);
    void Visualise(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in);
    void SetNumEigenvector(int e);
    void SetIteration(int i);
    void SetLambda(double l);
    void SetNoise(double n);
    
private:
    
    igl::opengl::glfw::Viewer& viewer;
    
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXd C;
    
    int iteration;
    double lambda;
    double noise;
    int eigenvector;
};
