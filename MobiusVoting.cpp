#include <igl/opengl/glfw/Viewer.h>
#include <igl/gaussian_curvature.h>
#include <igl/octree.h>
#include <igl/knn.h>
#include <iostream>
#include <ostream>
#include <algorithm>

#ifndef TRIANGLEMESHDS_HEADER
   #define TRIANGLEMESHDS_HEADER
   #include "TriangleMeshDS.cpp"
#endif

using namespace Eigen;
using namespace std;


class MobiusVoting
{
public:
    MobiusVoting( MatrixXd iV1, MatrixXd iF1, MatrixXd iV2, MatrixXd iF2){
        V1 = iV1;
        F1 = iF1;
        V2 = iV2;
        F2 = iF2;
    };
    ~MobiusVoting();

    void computeCorrespondances(int votingIterations, int minimalSubsetSize){
        // sampling potential correspondances for each mesh
        MatrixXd S1 = sample_correspondances(V1, F1, 50);
        MatrixXd S2 = sample_correspondances(V2, F2, 50);

        //

        for (int i = 0; i < votingIterations; i++){
            // take 3 random points from potential correspondances for each mesh
            VectorXi Z1(3);
            VectorXi Z2(3);
            do{
                for (int j = 0; j < 3; j++){
                    Z1(j) = std::rand() % S1.rows();
                    Z2(j) = std::rand() % S2.rows();
                }
            } while(Z1(0) == Z1(1) || Z1(2) == Z1(1) || Z1(0) == Z1(2) || Z2(0) == Z2(1) || Z2(2) == Z2(1) || Z2(0) == Z2(2));

            // compute mobius transformations for both triplets



        }


    }
private:
    MatrixXd V1, F1, V2, F2;

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
 * 
 * */
    MatrixXi mutual_nearest_neighbor(MatrixXd &set1, MatrixXd &set2){
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

        MatrixXi neighbors_2;
        std::vector<std::vector<int > > O_PI_2;
        MatrixXi O_CH_2;
        MatrixXd O_CN_2;
        VectorXd O_W_2;
        igl::octree(set2,O_PI_2,O_CH_2,O_CN_2,O_W_2);
        igl::knn(set2, set1, n_neighbors, O_PI_2, O_CH_2, O_CN_2, O_W_2, neighbors_2);

        // building MVN matrix 
        MatrixXi MVN = MatrixXi::Zero(n, n_neighbors);
        for (int i = 0; i < n; i++){

            for (int j = 0; j < n_neighbors; j++){
                bool is_neighbor = false;
                int k = 0;
                while (!is_neighbor && k < n_neighbors){
                    if (neighbors_1(i, j) == neighbors_2(j, k)){
                        is_neighbor = true;
                        MVN(i,j) += j + k;
                    }
                    k++;
                }
                if (!is_neighbor) MVN(i,j) = 100;
            }
        }

        // finding mutual neighbors
        std::vector<Vector3i> candidates;
        for (int i = 0; i < n; i++){

            for (int j = 0; j < n_neighbors; j++){
                Vector3i elt;
                elt << i, j, MVN(i,j);
                candidates.push_back(elt);
            }
        }
        std::sort(candidates.begin(), candidates.end(),
            [](Vector3i a, Vector3i b) {
                return (a(2) < b(2));
            } 
        );
    }    

};

