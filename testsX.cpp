#include "MobiusVoting.cpp"

using namespace Eigen; 


static int testsX(){
    MatrixXd V_test; 
    MatrixXi F_test;
    igl::readOFF("../data/cat0.off", V_test, F_test);

    MatrixXd V_test_alt = V_test + 0.1 * MatrixXd::Ones(V_test.rows(), V_test.cols()); 
    MobiusVoting votetest(V_test, F_test, V_test_alt, F_test);
    votetest.testMNN(V_test, V_test_alt);
}