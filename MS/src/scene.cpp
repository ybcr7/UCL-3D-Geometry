#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/file_exists.h>
#include "scene.h"
#include "ms.h"

Scene::Scene(igl::opengl::glfw::Viewer& refViewer):viewer(refViewer){
    default_C << 1.0,1.0,0.0;
}

Scene::~Scene(){}

void Scene::Discretisation(int mode){
    
    switch (mode){
        case 0:
			C_curvature = MS::UniformMeanCurvature(V,F);
            break;
        case 1:
			C_curvature = MS::GaussianCurvature(V,F);
            break;
        case 2:
			C_curvature = MS::NonUniformMeanCurvature(V,F);
            break;
        default:
            std::cout << "ERROR: Undefined Discretisation Mode" << std::endl;
            break;
    }

	VisualiseCurvature();

}

void Scene::Reconstruction() {
	ResetColor();
    Eigen::MatrixXd V_out(V.rows(),V.cols());
    V_out = MS::Reconstruction(V,F,eigenvector);
    Visualise(V_out, F);
}

void Scene::Smoothing(int mode) {
	ResetColor();
	V_smoothed = V_unsmoothed;
    Eigen::VectorXd C_out;
    switch (mode){
        case 0:
			V_smoothed = MS::ExplicitSmoothing(V_smoothed,F,lambda,iteration);
			C_curvature = MS::NonUniformMeanCurvature(V_smoothed, F);
            Visualise(V_smoothed,F);
            break;
        case 1:
			V_smoothed = MS::ImplicitSmoothing(V_smoothed,F,lambda,iteration);
			C_curvature = MS::NonUniformMeanCurvature(V_smoothed, F);
            Visualise(V_smoothed,F);
            break;
        default:
            std::cout << "ERROR: Undefined Smoothing Mode" << std::endl;
            break;
    }
	VisualiseCurvature();
}

void Scene::Initialise(std::string filename){

    std::vector<const char *> FILE_PATH_LIST{
            "../data/",
            "../../data/",
            "../../../data/"};
    bool file_found = false;

    for (const auto & FILE_PATH : FILE_PATH_LIST)
    {
        if ( igl::file_exists(FILE_PATH + filename) )
        {
            if (igl::readOFF(FILE_PATH+filename, V, F)) {
                file_found = true;
				V_unsmoothed = V;
				V_smoothed = V;
				C_curvature.resize(V.rows());
				C_curvature.setZero();
                break;
            }
        }
    }

    if (!file_found) {
        std::cout << "ERROR: "<< filename << " not found" << std::endl;
        exit(1);
    }

	ResetColor();
    Visualise(V,F);
}

void Scene::VisualiseComparison(int mode) {

	switch (mode) {
	case 0:
		Visualise(V_unsmoothed, F);
		break;
	case 1:
		Visualise(V_smoothed, F);
		break;
	default:
		Visualise(V, F);
		break;
	}
}

void Scene::AddNoise(){
	V_unsmoothed = MS::AddNoise(V, noise);
	V_smoothed = V_unsmoothed;
    Visualise(V_unsmoothed, F);
}

void Scene::SetNumEigenvector(int e) {
    if (e < 1) {
        eigenvector = 1;
    } else {
        eigenvector = e;
    }
}

void Scene::SetIteration(int i) {
    if (i < 1) {
        iteration = 1;
    } else {
        iteration = i;
    }
}

void Scene::SetLambda(double l) {
    if (l < 0) {
        lambda = 0;
    } else {
        lambda = l;
    }
}

void Scene::SetNoise(double n) {
    if (n < 0) {
        noise = 0;
    } else {
        noise = n;
    }
}

void Scene::SetCurvatureDisplayScale(double s) {
	if (s < 0) {
		curvature_display_scale = 0;
	}
	else {
		curvature_display_scale = s;
	}
}

void Scene::ResetColor() {
	C.resize(V.rows(), V.cols());
	for (int i = 0; i < V.rows(); i++) {
		C.row(i) = default_C;
	}
}

void Scene::Visualise(Eigen::MatrixXd V_in, Eigen::MatrixXi F_in){
    viewer.data().clear();
    viewer.data().show_overlay_depth = 1;
    viewer.data().set_mesh(V_in, F_in);
    viewer.data().set_colors(C);
    viewer.core.align_camera_center(V_in, F_in);
}

void Scene::VisualiseCurvature() {
	Eigen::VectorXd C_curvature_scaled;
	C_curvature_scaled = curvature_display_scale * C_curvature.array() / (C_curvature.maxCoeff() - C_curvature.minCoeff());
	igl::parula(C_curvature_scaled, false, C);
	Visualise(V_smoothed, F);

}

