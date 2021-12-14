#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/writeOBJ.h>
#include <igl/gaussian_curvature.h>
#include <igl/octree.h>
#include <igl/knn.h>
#include <igl/edges.h>
#include <iostream>
#include <ostream>

#include "MidEdgeDS.cpp"
#include "HarmonicSolver.cpp"
//#include "MobiusVoting.cpp"
#include "testsX.cpp"

using namespace Eigen; // to use the classes provided by Eigen library

MatrixXd V1; // matrix storing vertex coordinates of the input mesh (n rows, 3 columns)
MatrixXi F1; // incidence relations between faces and edges (f columns)

/**
 * This function is called every time a keyboard button is pressed
 * */
bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier) {
  std::cout << "pressed Key: " << key << " " << (unsigned int)key << std::endl;

  if (key == '1') {
    std::cout << "saving to OBJ format" << std::endl;
    igl::writeOBJ("../data/converted_mesh.obj", V1, F1);
  }

  if (key == '2')
  {
    // TO BE COMPLETED

    V1 = V1 / 5; 

    viewer.data(0).clear(); // Clear should be called before drawing the mesh
    viewer.data(0).set_mesh(V1, F1); // update the mesh (both coordinates and faces)
    //viewer.core().align_camera_center(V1, F1);
  }

  return false;
}

void draw_bounding_box(igl::opengl::glfw::Viewer &viewer, const MatrixXd &V)
{
  // compute the corners of the bounding box
  Vector3d m = V.colwise().minCoeff();
  Vector3d M = V.colwise().maxCoeff();

  MatrixXd V_box(8, 3);  // Corners of the bounding box
  MatrixXi E_box(12, 2); // edges of the bounding box

  V_box << m(0), m(1), m(2),
      M(0), m(1), m(2),
      M(0), M(1), m(2),
      m(0), M(1), m(2),
      m(0), m(1), M(2),
      M(0), m(1), M(2),
      M(0), M(1), M(2),
      m(0), M(1), M(2);

  E_box << 0, 1,
      1, 2,
      2, 3,
      3, 0,
      4, 5,
      5, 6,
      6, 7,
      7, 4,
      0, 4,
      1, 5,
      2, 6,
      7, 3;

  viewer.append_mesh();
  viewer.data(1).add_points(V_box, Eigen::RowVector3d(1, 0, 0));

  for (unsigned i = 0; i < E_box.rows(); ++i) // Plot the edges of the bounding box
    viewer.data().add_edges(
        V_box.row(E_box(i, 0)),
        V_box.row(E_box(i, 1)),
        Eigen::RowVector3d(1, 0, 0));
}

void draw_dots(igl::opengl::glfw::Viewer &viewer, const MatrixXd &V)
{

  
  viewer.append_mesh();
  viewer.data(1).add_points(V, RowVectorXd::Random(3));
  //  viewer.data(1).add_points(V, Eigen::RowVector3d(1, 0, 0));
}
/**
 * This function gets the potential correspondances vertices from the mesh
 * */
MatrixXd sample_correspondances(const MatrixXd &V, const MatrixXi &F, int N){
  VectorXd K;

  igl::gaussian_curvature(V,F,K);

  // getting the local maxima of of gaussian curvature

  int n = V.rows();
  int n_neighbors = 50;
  MatrixXi neighbors;
  std::vector<int> localMaxima;

  std::vector<std::vector<int > > O_PI;
  MatrixXi O_CH;
  MatrixXd O_CN;
  VectorXd O_W;
  igl::octree(V,O_PI,O_CH,O_CN,O_W);
  igl::knn(V, n_neighbors, O_PI, O_CH, O_CN, O_W, neighbors);

  VectorXd localK(n_neighbors);
  for(int l = 0; l < n; l++){
    for (int m = 0; m < n_neighbors; m++){
      localK(m) = K(neighbors(l,m));
    }
    if (K(l) >= localK.maxCoeff()){
      localMaxima.push_back(l);
    }
  } 

  MatrixXd sortedV(localMaxima.size(), 3);
  for (int k = 0; k < localMaxima.size(); k++){
    sortedV.row(k) = V.row(localMaxima[k]);
  }
  return sortedV;
}
// showcase function
void showcase(int feature){
  igl::opengl::glfw::Viewer viewer; // create the 3d viewer

  switch (feature) {
    case 1 : { // showcasing mid-edge mesh construction
      igl::readOFF("../data/cat0.off", V1, F1);
      //  print the number of mesh elements
      std::cout << "Vertices: " << V1.rows() << std::endl;
      std::cout << "Faces:    " << F1.rows() << std::endl;

      TriangleMeshDS *mesh = new TriangleMeshDS(V1, F1);
      MatrixXd angles;
      angles = mesh->compute_angles();
      
      MidEdgeDS *midedgemesh = new MidEdgeDS(*mesh);

      V1 = midedgemesh->getVert();
      F1 = midedgemesh->getFaces();

      viewer.data().set_mesh(V1, F1);
    }

    case 2 : { // showcasing selection of potential correspondance points
      igl::readOFF("../data/cat0.off", V1, F1);
      //  print the number of mesh elements
      std::cout << "Vertices: " << V1.rows() << std::endl;
      std::cout << "Faces:    " << F1.rows() << std::endl;

      MobiusVoting voter;// = new MobiusVoting;
      MatrixXd sampledPoints = voter.sample_correspondances(V1, F1, 50);
      viewer.data().set_mesh(V1, F1);
      draw_dots(viewer, sampledPoints);

    }
    case 3 : {// showcasing mid-edge flattening
      igl::readOFF("../data/cat0.off", V1, F1);
      //  print the number of mesh elements
      std::cout << "Vertices: " << V1.rows() << std::endl;
      std::cout << "Faces:    " << F1.rows() << std::endl;
    }
    case 4 : { // showcasing mutual nearest neighbors
      MatrixXd V_test, V_test_alt; 
      MatrixXi F_test, F_test_alt;
      igl::readOFF("../data/star.off", V_test, F_test);
      igl::readOFF("../data/star_warp.off", V_test_alt, F_test_alt);

      //MatrixXd V_test_alt = V_test + 0.1 * MatrixXd::Ones(V_test.rows(), V_test.cols()); 
      MobiusVoting votetest(V_test, F_test, V_test_alt, F_test_alt);
      votetest.testMNN(V_test, V_test_alt);

      viewer.data().set_mesh(V_test, F_test);
      
    }
  }
  viewer.launch();
}
// ------------ main program ----------------
int main(int argc, char *argv[])
{
 //testsX();

  igl::readOFF("../data/cat1.off", V1, F1); // Load an input mesh in OFF format
  //igl::readOFF("/Users/victor/Documents/ENSTA/inf574/projet/dataToOFF/converted_data/cat0.off", V1, F1); // Load an input mesh in OFF format

  //  print the number of mesh elements
  std::cout << "Vertices: " << V1.rows() << std::endl;
  std::cout << "Faces:    " << F1.rows() << std::endl;

  TriangleMeshDS *mesh = new TriangleMeshDS(V1, F1);
  MatrixXd angles;
  angles = mesh->compute_angles();

  // std::cout<<angles<<"\n";

  //int cutoff_face = mesh->get_cut_face();
  int cutoff_face = 1157;
  std::cout<<"CUTOFF FAAAACE !!! "<<cutoff_face<<"\n";
  
  MidEdgeDS *midedgemesh = new MidEdgeDS(*mesh);

  std::cout << "Vertices_midEdge: " << midedgemesh->getVert().rows() << std::endl;
  std::cout << "Faces_midEdge:    " << midedgemesh->getFaces().rows() << std::endl;
  std::cout<<"here0\n";

  MatrixXi faces = midedgemesh->getFaces();
  std::cout<<"here0\n";
  MatrixXi edges;
  edges = mesh->getEdges();


  std::cout<<"here0\n";

  HarmonicSolver *HS = new HarmonicSolver(V1.rows(), F1.rows(), midedgemesh->getVert().rows(), midedgemesh->getFaces().rows(), cutoff_face, angles);

  std::cout<<"here3\n";

  VectorXcd complex_flattening;
  complex_flattening = HS->get_complex_flattening(faces, F1, edges);

  
  std::cout<<"here9\n";
  MatrixXd vflatten(complex_flattening.rows(), 2);
  for (int i = 0; i < complex_flattening.rows(); i++){
    vflatten(i, 0) = complex_flattening(i).real();
    vflatten(i, 1) = complex_flattening(i).imag();
  }

  V1 = vflatten;
  F1 = midedgemesh->getFaces();

  igl::opengl::glfw::Viewer viewer; // create the 3d viewer
  viewer.callback_key_down = &key_down; // for dealing with keyboard events
  viewer.data().set_mesh(V1, F1); // load a face-based representation of the input 3d shape
  // MatrixXd sampledPoints = sample_correspondances(V1, F1, 20);
  // draw_dots(viewer, sampledPoints); // draw the boundaing box (red edges and vertices)

  viewer.launch(); // run the editor
}
