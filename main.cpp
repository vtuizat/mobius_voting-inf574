#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/writeOBJ.h>
#include <igl/gaussian_curvature.h>
#include <igl/octree.h>
#include <igl/knn.h>
#include <igl/edges.h>
#include <iostream>
#include <ostream>

//#include "MidEdgeDS.cpp"
//#include "HarmonicSolver.cpp"
//#include "MobiusVoting.cpp"
#include "testsX.cpp"

using namespace Eigen; // to use the classes provided by Eigen library

MatrixXd V1; // matrix storing vertex coordinates of the input mesh (n rows, 3 columns)
MatrixXi F1; // incidence relations between faces and edges (f columns)


void draw_dots(igl::opengl::glfw::Viewer &viewer, const MatrixXd &V)
{

  
  viewer.append_mesh();
  viewer.data(1).add_points(V, RowVectorXd::Random(3));
  //  viewer.data(1).add_points(V, Eigen::RowVector3d(1, 0, 0));
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
      break;
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
      break;

    }
    case 3 : {// showcasing mid-edge flattening
      igl::readOFF("../data/cat0.off", V1, F1);
      //  print the number of mesh elements
      std::cout << "Vertices: " << V1.rows() << std::endl;
      std::cout << "Faces:    " << F1.rows() << std::endl;
      break;
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
      break;
      
    }
    case 5 : { // showcasing full algo
      MatrixXd V_cat0, V_cat1; 
      MatrixXi F_cat0, F_cat1;
      igl::readOFF("../data/cat0.off", V_cat0, F_cat0);
      igl::readOFF("../data/cat1.off", V_cat1, F_cat1);

      MobiusVoting voter;
      MatrixXd sampledPoints1 = voter.sample_correspondances(V_cat0, F_cat0, 50);
      MatrixXd sampledPoints2 = voter.sample_correspondances(V_cat1, F_cat1, 50);

      VectorXcd mappedSampledPoints1 = voter.computeMapping( V_cat0, F_cat0, 1157);
      VectorXcd mappedSampledPoints2 = voter.sample_correspondances(V_cat1, F_cat1, 1157);

      MatrixXd correspondances = voter.computeCorrespondanceMatrix(sampledPoints1, sampledPoints2, mappedSampledPoints1, mappedSampledPoints2, 100, 10);

      std::cout << "\ncorrespondances : \n" << correspondances;

      break;
    }
  }
  viewer.launch();
}
// ------------ main program ----------------
int main(int argc, char *argv[])
{
 //testsX();
  showcase(1);

  igl::readOFF("../data/star.off", V1, F1); // Load an input mesh in OFF format
  //igl::readOFF("/Users/victor/Documents/ENSTA/inf574/projet/dataToOFF/converted_data/cat0.off", V1, F1); // Load an input mesh in OFF format

  //  print the number of mesh elements
  std::cout << "Vertices: " << V1.rows() << std::endl;
  std::cout << "Faces:    " << F1.rows() << std::endl;
  

  igl::opengl::glfw::Viewer viewer; // create the 3d viewer
  viewer.data().set_mesh(V1, F1); // load a face-based representation of the input 3d shape
  // MatrixXd sampledPoints = sample_correspondances(V1, F1, 20);
  // draw_dots(viewer, sampledPoints); // draw the boundaing box (red edges and vertices)

  viewer.launch(); // run the editor
}
