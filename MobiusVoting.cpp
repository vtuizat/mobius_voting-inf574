#include <igl/opengl/glfw/Viewer.h>
#include <igl/gaussian_curvature.h>
#include <igl/octree.h>
#include <igl/knn.h>
#include <iostream>
#include <ostream>
#include <algorithm>
#include <complex>
#include <cmath>
#include <iostream>




#ifndef TRIANGLEMESHDS_HEADER
   #define TRIANGLEMESHDS_HEADER
   #include "TriangleMeshDS.cpp"
#endif

using namespace Eigen;
using namespace std;


class MobiusVoting
{
public:
    MobiusVoting( MatrixXd iV1, MatrixXi iF1, MatrixXd iV2, MatrixXi iF2){
        V1 = iV1;
        F1 = iF1;
        V2 = iV2;
        F2 = iF2;
    };
    //~MobiusVoting();

/**
 * This function samples the potential correspondances vertices from the mesh
 * V,F : mesh data
 * N number of neighbors used to compute local Gauss curvature maxima
 * returns matrix of coordinates of potential correspondances vertices 
 * */
    MatrixXd sample_correspondances(const MatrixXd &V, const MatrixXi &F, int N){
        VectorXd K;

        igl::gaussian_curvature(V,F,K);

        // getting the local maxima of of gaussian curvature

        int n = V.rows();
        int n_neighbors = N;
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
/**
 * This function computes the correspondance matrix between S1 and S2 with votingIterations votes
 * */
    MatrixXd computeCorrespondanceMatrix(MatrixXd S1, MatrixXd S2, int votingIterations, int minimalSubsetSize){

        // sampling potential correspondances for each mesh
        //MatrixXd S1 = sample_correspondances(V1, F1, 50);
        //MatrixXd S2 = sample_correspondances(V2, F2, 50);

        //compute mapping

        VectorXcd mappedS1, mappedS2;

        MatrixXd C(S1.rows(), S2.rows()); // correspondence matrix
        double epsilon = 0.001;

        for (int i = 0; i < votingIterations; i++){
            // pick 3 random points from potential correspondances for each mesh
            VectorXi Z1(3);
            VectorXi Z2(3);
            do{
                for (int j = 0; j < 3; j++){
                    Z1(j) = std::rand() % S1.rows();
                    Z2(j) = std::rand() % S2.rows();
                }
            } while(Z1(0) == Z1(1) || Z1(2) == Z1(1) || Z1(0) == Z1(2) || Z2(0) == Z2(1) || Z2(2) == Z2(1) || Z2(0) == Z2(2));

            // compute mobius transformations for both triplets
            Vector3cd origin1, origin2, target;
            origin1 << mappedS1(Z1(0)), mappedS1(Z1(1)), mappedS1(Z1(2));
            origin2 << mappedS2(Z2(0)), mappedS2(Z2(1)), mappedS2(Z2(2));
            std::complex<double> comp_j( -1 / 2, sqrt(3) / 2);
            target << 1.0, comp_j, conj(comp_j);
            Matrix2cd mobius1 = mobius_interpolate(origin1, target);
            Matrix2cd mobius2 = mobius_interpolate(origin2, target);

            // apply mobius transformations to S1 and S2
            VectorXcd mobiusMappedS1 = mobius_apply(mappedS1, mobius1);
            VectorXcd mobiusMappedS2 = mobius_apply(mappedS2, mobius2);

            // find mutually nearest neighbors between mobiusMappedS2 and mobiusMappedS2
            MatrixXi MNN = mutual_nearest_neighbor(CtoR3(mobiusMappedS1),CtoR3(mobiusMappedS1));
            int valid_couples = 0;
            for (int k = 0; k < MNN.rows(); k++){
                if(MNN.row(k)(0) != -1) valid_couples++;
            }

            // voting process
            
            if (valid_couples > minimalSubsetSize){
                // deformation energy
                double E_c = computeDeformationEnergy(CtoR3(mobiusMappedS1),CtoR3(mobiusMappedS1),MNN);
                for (int k = 0; k < MNN.rows(); k++){
                    if(MNN.row(k)(0) != -1) {
                        C(k,MNN.row(k)(0)) += 1 / (epsilon + E_c);
                    }
                }
            }

        }

        return C;
    }
/**
 * This function takes a correspondance matrix and returns a vector of coupled correspondances <mesh1, mesh2>  
 * */
    vector<pair<int,int>> processCorrespondanceMatrix(MatrixXd C){
        int n_rows = C.rows();
        int n_cols =C.cols(); 
        vector<pair<int,int>> correspondances;

        pair<int,int> maxIndex;
        double maxC = C.maxCoeff(&maxIndex.first, &maxIndex.second); 
        if (maxC != 0) C /= maxC;

        unsigned i = 0;
        while (maxC != 0 && i < C.size()){
            correspondances.push_back(maxIndex);
            C.row(maxIndex.first) = VectorXd::Zero(n_rows);
            C.col(maxIndex.second) = VectorXd::Zero(n_cols);
            maxC = C.maxCoeff(&maxIndex.first, &maxIndex.second);
            i++;
        }

        return correspondances;
    }

// TESTS
    void testMNN(MatrixXd set1, MatrixXd set2){
        cout << "set 1 :\n" << set1 << "\nset 2 :\n" << set2 <<"\n";
        MatrixXi M = mutual_nearest_neighbor(set1,set2);
        cout << "MNN :\n";
        for (int i =0; i < M.rows(); i++){
            cout << i << " " << M(i, 0) << "\n";
        }
    } 

private:
    MatrixXd V1, V2;
    MatrixXi F1, F2;

/**
 * This function computes the mobius transformation that maps origin to target
 * origin : complex triplet
 * target : complex triplet
 * returns matrix with complex mobius transformation parameters (a b, c d) 
 * */
    Matrix2cd mobius_interpolate(Vector3cd origin, Vector3cd target){
        //
        Matrix2cd A, B;
        A << target(1) - target(2), target(0) * target(2) - target(0) * target(1),
            target(1) - target(0), target(0) * target(2) - target(2) * target(1);
        B << origin(1) - origin(2), origin(0) * origin(2) - origin(0) * origin(1),
            origin(1) - origin(0), origin(0) * origin(2) - origin(2) * origin(1);
        Matrix2cd mobius = A.inverse() * B;
        return mobius;
    }
/**
 * This function applies a mobius transformation to each point from a list of points 
 * and returns a list of transformed points
 * */
    VectorXcd mobius_apply(VectorXcd points, Matrix2cd mobius){
        VectorXcd transformed_points(points.size());

        for (int k = 0; k < points.size(); k++){
            transformed_points(k) = (mobius(0,0) * points(k) + mobius(0,1)) / (mobius(1,0) * points(k) + mobius(1,1));
        }
        return transformed_points;
    }

/**
 * returns the correspondances form set2 for each element of set1, if no correspondqnce the value is set to -1
 * */
    MatrixXi mutual_nearest_neighbor(MatrixXd set1, MatrixXd set2){
        int n = set1.rows();

        // knn from the other set for set1 and set2
        int n_neighbors = 5;
        MatrixXi neighbors_1;
        std::vector<std::vector<int > > O_PI_1;
        MatrixXi O_CH_1;
        MatrixXd O_CN_1;
        VectorXd O_W_1;
        igl::octree(set1,O_PI_1,O_CH_1,O_CN_1,O_W_1);
        igl::knn(set1, set2, n_neighbors, O_PI_1, O_CH_1, O_CN_1, O_W_1, neighbors_1);
        cout << neighbors_1 << "\n";
        MatrixXi neighbors_2;
        std::vector<std::vector<int > > O_PI_2;
        MatrixXi O_CH_2;
        MatrixXd O_CN_2;
        VectorXd O_W_2;
        igl::octree(set2,O_PI_2,O_CH_2,O_CN_2,O_W_2);
        igl::knn(set2, set1, n_neighbors, O_PI_2, O_CH_2, O_CN_2, O_W_2, neighbors_2);
        cout << neighbors_2 << "\n";
        // building MVN matrix 
        MatrixXi MVN = MatrixXi::Zero(n, n_neighbors);
        for (int i = 0; i < n; i++){

            for (int j = 0; j < n_neighbors; j++){
                /* bool is_neighbor = false;
                int k = 0;
                while (!is_neighbor && k < set2.rows()){
                    if (neighbors_1(i, j) == neighbors_2(k, j)){
                        is_neighbor = true;
                        MVN(i,j) += j + k;
                    }
                    k++;
                }
                if (!is_neighbor) MVN(i,j) = 100; */
                MVN(i,j) = mvn(i, neighbors_1(i,j), neighbors_1, neighbors_2);
            }
        }

        cout << "\n MVN mat : \n" << MVN;

        // finding mutual neighbors
        std::vector<Vector3i> candidates, c_copy;
        
        for (int i = 0; i < n; i++){

            for (int j = 0; j < n_neighbors; j++){
                Vector3i elt;
                elt << i, neighbors_1(i,j), MVN(i,j);
                cout << "\n elt : " << elt.transpose() << "\n";
                candidates.push_back(elt);
            }
        }
        c_copy = candidates;
        cout << "\n candidates : \n" ;
        for (int o = 0; o<candidates.size(); o++) cout << candidates.at(o)(0) << " " << candidates.at(o)(1) << " " << candidates.at(o)(2) << "\n" ;
        std::sort(candidates.begin(), candidates.end(),
            [](Vector3i a, Vector3i b) {
                return (a(2) < b(2));
            } 
        );
        /* bool is_equal = candidates == c_copy;
        cout << "\n sorted : \n";
        for (int o = 0; o<candidates.size(); o++) cout << candidates.at(o)(0) << " " << candidates.at(o)(1) << " " << candidates.at(o)(2) << "\n" ;
        cout << "equal" << is_equal; */

        //selecting the best pairs
        MatrixXi correspondances = -1 * MatrixXi::Ones(n, 2);
        int count = 0;
        int paired_count = 0;
        while (count < candidates.size() && paired_count < n){
            if(correspondances(candidates.at(count)(0),0) == -1 && candidates.at(count)(2) != 100){
                correspondances.row(candidates.at(count)(0)) << candidates.at(count)(1), candidates.at(count)(2);
                paired_count++;
            }
            else if(candidates.at(count)(2) == correspondances.row(candidates.at(count)(0))(1) && candidates.at(count)(2) != 100 &&
            (set1.row(candidates.at(count)(0)) - set2.row(candidates.at(count)(1))).norm() < (set1.row(candidates.at(count)(0)) - set2.row(correspondances.row(candidates.at(count)(0))(0) )).norm() ){
                correspondances.row(candidates.at(count)(0)) << candidates.at(count)(1), candidates.at(count)(2);
            }
            count++;
            //cout << "\n count = " << count << "\n paired_count = " << paired_count <<" coresp : \n" << correspondances <<"\n";
        }

        return correspondances.col(0);
    }
/**
 * computes the mutual neighborhood value for i, j
 * */
    int mvn(int i, int j, MatrixXi neighbors_1,MatrixXi neighbors_2){
        int value = 100;
        for (int x = 0; x < neighbors_1.cols(); x++){
            if(neighbors_1(i, x) == j) value = x;
        }
        for (int x = 0; x < neighbors_2.cols(); x++){
            if(neighbors_2(j, x) == i && value != 100) value += x;
        }
        return value;
    }

/**
 * transforms a list of complex points into a list of points in R3 (z=0)
 * */
    MatrixXd CtoR3(VectorXcd points){
        MatrixXd M(points.size(), 3);
        for (int i = 0; i < points.size(); i++){
            M.row(i) << points(i).real(), points(i).imag(), 0.0;
        }
        return M;
    }
    
/**
 * computes deformation energy for a given set of correspondances
 *  
 * */
    double computeDeformationEnergy(MatrixXd set1, MatrixXd set2, MatrixXi MNN){
        double E_c = 0.0;
        for (int k = 0; k < MNN.rows(); k++){
            if(MNN.row(k)(0) != -1) {
                E_c += (set1.row(k)-set2.row(MNN.row(k)(0))).norm();
            }
        }
        return E_c;
    }



};

