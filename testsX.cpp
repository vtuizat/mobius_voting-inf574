#include "MobiusVoting.cpp"

using namespace Eigen; 


static int testsX(){
    MatrixXd V_test, V_test_alt; 
    MatrixXi F_test, F_test_alt;
    igl::readOFF("../data/star.off", V_test, F_test);
    igl::readOFF("../data/star_warp.off", V_test_alt, F_test_alt);

    //MatrixXd V_test_alt = V_test + 0.1 * MatrixXd::Ones(V_test.rows(), V_test.cols()); 
    MobiusVoting votetest(V_test, F_test, V_test_alt, F_test_alt);
    votetest.testMNN(V_test, V_test_alt);
    return 0;
}