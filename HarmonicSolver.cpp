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
    HarmonicSolver(int nV, int nF, int nV_ME, int nF_ME, int cutoff_face_in, MatrixXd  &angles_in){
        
        nVertices = nV;
        nFaces = nF;
        nVertices_ME = nV_ME;
        nFaces_ME = nF_ME;

        cutoff_face = cutoff_face_in;

        angles = &angles_in;
        //harmonic_weights = new MatrixXd(nVertices_ME, 2);
        
        std::cout<<"here1\n";
        harmonic_weights = new MatrixXd(nVertices_ME, 2);
        std::cout<<"here1\n";
        (*harmonic_weights) = MatrixXd::Constant(nVertices_ME, 2, -1.0);
        
        std::cout<<"here2\n";

        //cut_off_face = min_geodesic_distance();

    }

    VectorXcd get_complex_flattening(MatrixXi &faces, MatrixXi &F, MatrixXi &E){
        VectorXcd res(nVertices_ME);
        std::cout<<"here4\n";
        std::complex<double> icomplex(0.0, 1.0);
        compute_harmonic_weights(faces, F, E);
        std::cout<<"here7\n";
        for (int i = 0; i<nVertices_ME; i++){
            res(i) = (*harmonic_weights)(i, 0) +icomplex*(*harmonic_weights)(i, 1);
        }
        std::cout<<"here8\n";
        return res;
    }


    
private:

        void compute_harmonic_weights(MatrixXi &faces, MatrixXi &F, MatrixXi &E){

        int v1 = F(cutoff_face, 0);
        int v2 = F(cutoff_face, 1);
        
        std::cout<<"here5\n";

        MatrixXd L = compute_matrix(F);

        VectorXd rhs = VectorXd::Zero(nVertices);

        VectorXd zeros = VectorXd::Zero(nVertices);
        L.row(v1) = zeros;
        L.row(v2) = zeros;
        L(v1, v1) = 1.0;
        L(v2, v2) = 1.0;

        rhs(v1) = 1.0;
        rhs(v2) = -1.0;
        
        std::cout<<"here6\n";

        // FullPivLU<MatrixXd> solver;
        // solver.compute(L);
        

        SparseMatrix<double> A(nVertices, nVertices);
        A = L.sparseView();
        // A.makeCompressed();

        SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> solver;

        solver.compute(A);
        if(solver.info()!=Success) {
            std::cout<<"Decomposition failed !\n";
            return;
        }

        // solver.analyzePattern(A);
        // if(solver.info()!=Success) {
        //     std::cout<<"Decomposition failed !\n";
        //     return;
        // }

        // solver.factorize(A);
        // if(solver.info()!=Success) {
        //     std::cout<<"Decomposition failed !\n";
        //     return;
        // }

        VectorXd res(nVertices);
        res = solver.solve(rhs);
        if(solver.info()!=Success) {
            std::cout<<"Solving failed !\n";
            return;
        }
        std::cout<<"here6\n";
        std::cout<<res<<"\n";

        int i, j;
        for (int edge; edge < nVertices_ME; edge++){
            i = E(edge, 0);
            j = E(edge, 1);
            (*harmonic_weights)(edge, 0) = (res(i) + res(j))/2;

        }
        std::cout<<"here6\n";  

        compute_conjugates(faces, F, E, res);

    }

    void compute_conjugates(MatrixXi &faces, MatrixXi &F, MatrixXi &E, VectorXd &u){

        (*harmonic_weights)(0, 1) = 0;
        
        int n_edges = 0;
        Vector2i edge1, edge2;
        int v1, v2, v3;
        int face;
        int last, commun;

        while(n_edges < nVertices_ME-1){
            std::cout<<"Number of conjugates : "<<n_edges<<"\n";
            for (int face_ME = 0; face_ME < nFaces_ME; face_ME++) {

                v1 = faces(face_ME, 0);
                v2 = faces(face_ME, 1);
                v3 = faces(face_ME, 2);

                if ((*harmonic_weights)(v1, 1) == -1 && (*harmonic_weights)(v2, 1) != -1){
                    edge1 = E.row(v1);
                    edge2 = E.row(v2);
                    last = 0;
                    if (edge2(last) == edge1(0))  {last = 1; commun = 0;}
                    else if (edge2(last) == edge1(1))  {last = 1; commun = 1;}
                    face = getFaceIndexFromVert(edge1(0), edge1(1), edge2(last), F);

                    (*harmonic_weights)(v1, 1) = (*harmonic_weights)(v2, 1)
                                                 + ((u(edge1(!commun))-u(edge1(commun)))*cot(getAngleFromFaceAndVertex(face, edge2(last), F)) 
                                                 + (u(edge2(last))-u(edge2(!last)))*cot(getAngleFromFaceAndVertex(face, edge1(!commun), F)))/2;
                    n_edges ++;

                }

                else if ((*harmonic_weights)(v2, 1) == -1 && (*harmonic_weights)(v1, 1) != -1){
                    edge1 = E.row(v2);
                    edge2 = E.row(v1);
                    last = 0;
                    if (edge2(last) == edge1(0))  {last = 1; commun = 0;}
                    else if (edge2(last) == edge1(1))  {last = 1; commun = 1;}
                    face = getFaceIndexFromVert(edge1(0), edge1(1), edge2(last), F);

                    (*harmonic_weights)(v2, 1) = (*harmonic_weights)(v1, 1)
                                                 + ((u(edge1(!commun))-u(edge1(commun)))*cot(getAngleFromFaceAndVertex(face, edge2(last), F)) 
                                                 + (u(edge2(last))-u(edge2(!last)))*cot(getAngleFromFaceAndVertex(face, edge1(!commun), F)))/2;
                    n_edges ++;

                }

                else if ((*harmonic_weights)(v2, 1) == -1 && (*harmonic_weights)(v3, 1) != -1){
                    edge1 = E.row(v2);
                    edge2 = E.row(v3);
                    last = 0;
                    if (edge2(last) == edge1(0))  {last = 1; commun = 0;}
                    else if (edge2(last) == edge1(1))  {last = 1; commun = 1;}
                    face = getFaceIndexFromVert(edge1(0), edge1(1), edge2(last), F);

                    (*harmonic_weights)(v2, 1) = (*harmonic_weights)(v3, 1)
                                                 + ((u(edge1(!commun))-u(edge1(commun)))*cot(getAngleFromFaceAndVertex(face, edge2(last), F)) 
                                                 + (u(edge2(last))-u(edge2(!last)))*cot(getAngleFromFaceAndVertex(face, edge1(!commun), F)))/2;
                    
                    n_edges ++;
                    
                }

                else if ((*harmonic_weights)(v3, 1) == -1 && (*harmonic_weights)(v2, 1) != -1){
                    edge1 = E.row(v3);
                    edge2 = E.row(v2);
                    last = 0;
                    if (edge2(last) == edge1(0))  {last = 1; commun = 0;}
                    else if (edge2(last) == edge1(1))  {last = 1; commun = 1;}
                    face = getFaceIndexFromVert(edge1(0), edge1(1), edge2(last), F);

                    (*harmonic_weights)(v3, 1) = (*harmonic_weights)(v2, 1)
                                                 + ((u(edge1(!commun))-u(edge1(commun)))*cot(getAngleFromFaceAndVertex(face, edge2(last), F)) 
                                                 + (u(edge2(last))-u(edge2(!last)))*cot(getAngleFromFaceAndVertex(face, edge1(!commun), F)))/2;
                    n_edges ++;
                    
                }

                else if ((*harmonic_weights)(v3, 1) == -1 && (*harmonic_weights)(v1, 1) != -1){
                    edge1 = E.row(v3);
                    edge2 = E.row(v1);
                    last = 0;
                    if (edge2(last) == edge1(0))  {last = 1; commun = 0;}
                    else if (edge2(last) == edge1(1))  {last = 1; commun = 1;}
                    face = getFaceIndexFromVert(edge1(0), edge1(1), edge2(last), F);

                    (*harmonic_weights)(v3, 1) = (*harmonic_weights)(v1, 1)
                                                 + ((u(edge1(!commun))-u(edge1(commun)))*cot(getAngleFromFaceAndVertex(face, edge2(last), F)) 
                                                 + (u(edge2(last))-u(edge2(!last)))*cot(getAngleFromFaceAndVertex(face, edge1(!commun), F)))/2;
                    n_edges ++;
                    
                }

                else if ((*harmonic_weights)(v1, 1) == -1 && (*harmonic_weights)(v3, 1) != -1){
                    edge1 = E.row(v1);
                    edge2 = E.row(v3);
                    last = 0;
                    if (edge2(last) == edge1(0))  {last = 1; commun = 0;}
                    else if (edge2(last) == edge1(1))  {last = 1; commun = 1;}
                    face = getFaceIndexFromVert(edge1(0), edge1(1), edge2(last), F);

                    (*harmonic_weights)(v1, 1) = (*harmonic_weights)(v3, 1)
                                                 + ((u(edge1(!commun))-u(edge1(commun)))*cot(getAngleFromFaceAndVertex(face, edge2(last), F)) 
                                                 + (u(edge2(last))-u(edge2(!last)))*cot(getAngleFromFaceAndVertex(face, edge1(!commun), F)))/2;
                    n_edges ++;
                    
                }
            } 
        }

    }

    // MatrixXd compute_matrix_triplets(MatrixXi &F){
        
    //     //MatrixXd L = MatrixXd::Zero(nVertices, nVertices);
    //     std::vector<Triplet> L;

    //     for (int i = 0; i < nVertices; i++){

    //         float sum = 0;

    //         MatrixXd adjacency_angles = compute_adjacency(i, F);

    //         for (int j = 0; j < adjacency_angles.rows(); j++){

    //             int adj_vert = int(adjacency_angles(j, 0));
    //             if(adjacency_angles(j, 2) == -1.0) {
    //                 std::cout<<"Edge at "<<i<<"-"<<adj_vert<<"\n";
    //                 L.push_back(Triplet(i, adj_vert, cot(adjacency_angles(j, 1))/2.0));
    //             }
    //             else {
    //                 L.push_back(Triplet(i, adj_vert, (cot(adjacency_angles(j, 1))+cot(adjacency_angles(j, 2)))/2.0));
    //             }
    //             sum += L[-1].third;
    //         }

    //         L(i, i) = -sum;
    //     }

    //     return L;

    // }

    MatrixXd compute_matrix(MatrixXi &F){
        
        MatrixXd L = MatrixXd::Zero(nVertices, nVertices);

        for (int i = 0; i < nVertices; i++){

            float sum = 0;

            MatrixXd adjacency_angles = compute_adjacency(i, F);

            for (int j = 0; j < adjacency_angles.rows(); j++){

                int adj_vert = int(adjacency_angles(j, 0));
                if(adjacency_angles(j, 2) == -1.0) {
                    std::cout<<"Edge at "<<i<<"-"<<adj_vert<<"\n";
                    L(i, adj_vert) = cot(adjacency_angles(j, 1))/2.0;
                    std::cout<<"Adj : "<<adjacency_angles(j, 1)<<"\n";
                    std::cout<<"MatrixValue : "<<L(i, adj_vert)<<"\n";
                }
                else {
                    L(i, adj_vert) = (cot(adjacency_angles(j, 1))+cot(adjacency_angles(j, 2)))/2.0;
                    std::cout<<"Adj : "<<adjacency_angles(j, 1)<<"\n";
                    std::cout<<"Adj : "<<adjacency_angles(j, 2)<<"\n";
                    std::cout<<"MatrixValue : "<<L(i, adj_vert)<<"\n";
                }
                sum += L(i, adj_vert);
            }

            L(i, i) = -sum;
            std::cout<<"MatrixDiagValue : "<<L(i, i)<<"\n";

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

	int getFaceIndexFromVert(int v1, int v2, int v3, MatrixXi &F){
		
		bool nf = false;
		int f = 0;
		bool v1_in_face, v2_in_face, v3_in_face;

		while (!nf && f < nFaces) {

			v1_in_face = (F.row(f).array() == v1).any();
			v2_in_face = (F.row(f).array() == v2).any();
			v3_in_face = (F.row(f).array() == v3).any();
			if (v1_in_face && v2_in_face && v3_in_face){
				nf = true;
                return f;
			}
			f++;
		}

		return f;
	}

    int getAngleFromFaceAndVertex(int face, int vertex, MatrixXi &F){
        int idx = 0;
        while (F(face, idx) != vertex) idx++;
        return (*angles)(face, idx);
    }

    double cot(double a1){
        return cos(a1)/sin(a1);
    }

    int nVertices, nFaces, nVertices_ME, nFaces_ME;
    int cutoff_face;
    MatrixXd *angles;
    MatrixXd *harmonic_weights;

};