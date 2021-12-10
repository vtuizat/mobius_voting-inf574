#include <igl/opengl/glfw/Viewer.h>
#include <igl/vertex_triangle_adjacency.h>
#include <iostream>
#include <ostream>

#ifndef TRIANGLEMESHDS_HEADER
   #define TRIANGLEMESHDS_HEADER
   #include "TriangleMeshDS.cpp"
#endif

using namespace Eigen;
using namespace std;

class HarmonicSolver
{

public:
    HarmonicSolver(int nV, int nF, MatrixXd  &angles_in){
        
        nVertices = nV;
        nFaces = nF;
        angles = &angles_in;
        harmonic_weights = new MatrixXd(nVertices, 2);

        //cut_off_face = min_geodesic_distance();

    }

    void compute_harmonic_weights(MatrixXi &F){

        MatrixXd L = compute_matrix(F);

        //solve L;

        compute_conjugates(F);
        
    }
    
private:

    void compute_conjugates(MatrixXi &F){
        (*harmonic_weights)(F(0, 0)) = 0;
        std::vector<std::vector<int>> VF, VI;
        igl::vertex_triangle_adjacency(nVertices, F, VF, VI);
        for (int i = 0; i < nVertices; i++){
            for (int j = 0; j < VF[i].size(); j++){
                std::cout<<VF[i][j]<<", ";
            }
            std::cout<<"\n";
        }
    }

    MatrixXd compute_matrix(MatrixXi &F){
        
        MatrixXd L(nVertices, nVertices);

        for (int i = 0; i < nVertices; i++){

            float sum = 0;

            MatrixXd adjacency_angles = compute_adjacency(i, F);

            for (int j = 0; j < adjacency_angles.rows(); j++){

                int adj_vert = int(adjacency_angles(j, 0));
                if(adjacency_angles(j, 2) == -1.0) {
                    std::cout<<"Edge at "<<i<<"-"<<adj_vert<<"\n";
                    L(i, adj_vert) = cot(adjacency_angles(j, 1))/2.0;
                }
                else {
                    L(i, adj_vert) = (cot(adjacency_angles(j, 1))+cot(adjacency_angles(j, 2)))/2.0;
                }
                sum += L(i, adj_vert);
            }

            L(i, i) = -sum;
        }

        return L;

    }

    MatrixXd compute_adjacency(int v, MatrixXi &F){

        VectorXi adj;
        VectorXi adj_faces;

        get1Ring(v, F, adj, adj_faces);

        int nAdjVertices = adj.size();
        int nAdjFaces = adj_faces.size();


        MatrixXd support_angles = MatrixXd::Constant(nAdjVertices, 3, -1.0);
        VectorXd adj_d;
        adj_d = adj.cast<double>();
        support_angles.col(0) = adj_d;
        
        get_support_angles(v, F, adj, adj_faces, support_angles);
        
        return support_angles;

    }

   void get1Ring(int v, MatrixXi &F, VectorXi &adj, VectorXi &faces_idx){

        std::vector<int> adjacency;
        std::vector<int> adjacency_face;

        int nAdjFaces = 0;

        for (int face = 0; face < nFaces; face++){

            if ((F.row(face).array() == v).any()){

                adjacency_face.push_back(face);
                adjacency.push_back(F(face, 0));
                adjacency.push_back(F(face, 1));
                adjacency.push_back(F(face, 2));

                nAdjFaces++;
            }
        }

        std::sort(adjacency.begin(), adjacency.end());
        adjacency.erase(std::unique(adjacency.begin(), adjacency.end()), adjacency.end());
        adjacency.erase(std::remove(adjacency.begin(), adjacency.end(), v), adjacency.end());

        adj = Map<VectorXi, Unaligned>(adjacency.data(), adjacency.size());
        faces_idx = Map<VectorXi, Unaligned>(adjacency_face.data(), adjacency_face.size());
    }

    void get_support_angles(int v, MatrixXi &F, VectorXi adj, VectorXi adj_faces, MatrixXd &support_angles){

        for (int id = 0; id < adj_faces.size(); id++){

            int face = adj_faces(id);
            
            int vertex = 0;
            while (F(face, vertex) != v)vertex++;

            int opp_vert1 = F(face, (vertex+1)%3);
            int opp_vert2 = F(face, (vertex+2)%3);

            int idx1 = get_index(adj, opp_vert1);
            int idx2 = get_index(adj, opp_vert2);

            if (support_angles(idx1, 1) == -1.0) support_angles(idx1, 1) = (*angles)(face, (vertex+2)%3);
            else support_angles(idx1, 2) = (*angles)(face, (vertex+2)%3);

            if (support_angles(idx2, 1) == -1.0) support_angles(idx2, 1) = (*angles)(face, (vertex+1)%3);
            else support_angles(idx2, 2) = (*angles)(face, (vertex+1)%3);

        }

    }

    int get_index(VectorXi v, int value){

        bool test = true;
        int index = 0;
        while (v(index) != value)index++;

        return index;
    }

    double cot(double a1){
        return cos(a1)/sin(a1);
    }

    int nVertices, nFaces;
    int cut_off_face;
    MatrixXd *angles;
    MatrixXd *harmonic_weights;

};