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
    HarmonicSolver(int nV, int nF, int nV_ME, int nF_ME, int cutoff_face_in, MatrixXd &angles_in){
        
        nVertices = nV;
        nFaces = nF;
        nVertices_ME = nV_ME;
        nFaces_ME = nF_ME;

        cutoff_face = cutoff_face_in;

        angles = &angles_in;
        //harmonic_weights = new MatrixXd(nVertices_ME, 2);
        
        harmonic_weights = new MatrixXd(nVertices_ME, 2);
        (*harmonic_weights) = MatrixXd::Constant(nVertices_ME, 2, -1.0);

        //cut_off_face = min_geodesic_distance();

    }

    VectorXcd get_complex_flattening(MatrixXi &faces, MatrixXi &F, MatrixXi &E){
        VectorXcd res(nVertices_ME);
        const std::complex<double> icomplex(0.0, 1.0);
        compute_harmonic_weights(faces, F, E);
        for (int i = 0; i<nVertices_ME; i++){
            res(i) = (*harmonic_weights)(i, 0) +icomplex*(*harmonic_weights)(i, 1);
        }
        return res;
    }

private:

        void compute_harmonic_weights(MatrixXi &faces, MatrixXi &F, MatrixXi &E){

        int v1 = F(cutoff_face, 0);
        int v2 = F(cutoff_face, 1);

        MatrixXd L = compute_matrix(F);

        VectorXd rhs = VectorXd::Zero(nVertices);

        VectorXd zeros = VectorXd::Zero(nVertices);
        L.row(v1) = zeros;
        L.row(v2) = zeros;
        L(v1, v1) = 1.0;
        L(v2, v2) = 1.0;

        rhs(v1) = 1.0;
        rhs(v2) = -1.0;

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
        
        int i, j;
        for (int edge = 0; edge < nVertices_ME; edge++){
            i = E(edge, 0);
            j = E(edge, 1);
            (*harmonic_weights)(edge, 0) = (res(i) + res(j))/2.0;

        }

        compute_conjugates_v2(faces, F, E, res);


    }

    void compute_conjugates_v2(MatrixXi &faces, MatrixXi &F, MatrixXi &E, VectorXd &u){

        std::cout<<"Computing conjugate harmonic weights...\n";
        
        (*harmonic_weights)(0, 1) = 0;

        VectorXi adj_faces;

        std::vector<int> source_edges;
        std::vector<int> next_edges;

        source_edges.resize(1);
        source_edges[0] = 0;
        
        int nE = 0;
        int edge1, edge2, face, opp1, opp2;
        double angle1, angle2;
        while (nE < nVertices_ME-1){
            for (int se = 0; se < source_edges.size(); se++){
                edge1 = source_edges[se];
                VectorXi one_ring;
                get1Ring(edge1, faces, one_ring, adj_faces);
                for (int te = 0; te < one_ring.size(); te++){
                    edge2 = one_ring(te);
                    if ((*harmonic_weights)(edge2, 1) == -1.0){
                        compute_one_conjugate(edge1, edge2, F, E, u);
                        next_edges.push_back(edge2);
                        nE++;
                    }
                }
            }
            source_edges.clear();
            source_edges.resize(next_edges.size());
            source_edges = next_edges;
            next_edges.clear();
        }
    }


    void compute_one_conjugate(int source_edge, int target_edge, MatrixXi &F, MatrixXi &E, VectorXd &u){
        int opp1, opp2, sommet, face;
        double angle1, angle2;
        face = getFaceIndexFromEdge(E.row(source_edge), E.row(target_edge), F, opp1, opp2, sommet);
        angle1 = getAngleFromFaceAndVertex(face, opp1, F);
        angle2 = getAngleFromFaceAndVertex(face, opp2, F);
        (*harmonic_weights)(target_edge, 1) = (*harmonic_weights)(source_edge, 1)
                                              + ((u(opp2)-u(sommet))*cot(angle1) + (u(opp1)-u(sommet))*cot(angle2))/2.0;           

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

                    // std::cout<<getAngleFromFaceAndVertex(face, edge2(last), F)<<"\n";
                    // std::cout<<getAngleFromFaceAndVertex(face, edge2(!commun), F)<<"\n\n";

                    (*harmonic_weights)(v1, 1) = (*harmonic_weights)(v2, 1)
                                                 + ((u(edge1(!commun))-u(edge1(commun)))*cot(getAngleFromFaceAndVertex(face, edge2(last), F)) 
                                                 + (u(edge2(last))-u(edge2(!last)))*cot(getAngleFromFaceAndVertex(face, edge1(!commun), F)))/2.0;
                    n_edges ++;

                }

                else if ((*harmonic_weights)(v2, 1) == -1 && (*harmonic_weights)(v1, 1) != -1){
                    edge1 = E.row(v2);
                    edge2 = E.row(v1);
                    last = 0;
                    if (edge2(last) == edge1(0))  {last = 1; commun = 0;}
                    else if (edge2(last) == edge1(1))  {last = 1; commun = 1;}
                    face = getFaceIndexFromVert(edge1(0), edge1(1), edge2(last), F);
                    // std::cout<<getAngleFromFaceAndVertex(face, edge2(last), F)<<"\n";
                    // std::cout<<getAngleFromFaceAndVertex(face, edge2(!commun), F)<<"\n";

                    (*harmonic_weights)(v2, 1) = (*harmonic_weights)(v1, 1)
                                                 + ((u(edge1(!commun))-u(edge1(commun)))*cot(getAngleFromFaceAndVertex(face, edge2(last), F)) 
                                                 + (u(edge2(last))-u(edge2(!last)))*cot(getAngleFromFaceAndVertex(face, edge1(!commun), F)))/2.0;
                    n_edges ++;

                }

                else if ((*harmonic_weights)(v2, 1) == -1 && (*harmonic_weights)(v3, 1) != -1){
                    edge1 = E.row(v2);
                    edge2 = E.row(v3);
                    last = 0;
                    if (edge2(last) == edge1(0))  {last = 1; commun = 0;}
                    else if (edge2(last) == edge1(1))  {last = 1; commun = 1;}
                    face = getFaceIndexFromVert(edge1(0), edge1(1), edge2(last), F);
                    // std::cout<<getAngleFromFaceAndVertex(face, edge2(last), F)<<"\n";
                    // std::cout<<getAngleFromFaceAndVertex(face, edge2(!commun), F)<<"\n";

                    (*harmonic_weights)(v2, 1) = (*harmonic_weights)(v3, 1)
                                                 + ((u(edge1(!commun))-u(edge1(commun)))*cot(getAngleFromFaceAndVertex(face, edge2(last), F)) 
                                                 + (u(edge2(last))-u(edge2(!last)))*cot(getAngleFromFaceAndVertex(face, edge1(!commun), F)))/2.0;
                    
                    n_edges ++;
                    
                }

                else if ((*harmonic_weights)(v3, 1) == -1 && (*harmonic_weights)(v2, 1) != -1){
                    edge1 = E.row(v3);
                    edge2 = E.row(v2);
                    last = 0;
                    if (edge2(last) == edge1(0))  {last = 1; commun = 0;}
                    else if (edge2(last) == edge1(1))  {last = 1; commun = 1;}
                    face = getFaceIndexFromVert(edge1(0), edge1(1), edge2(last), F);
                    // std::cout<<getAngleFromFaceAndVertex(face, edge2(last), F)<<"\n";
                    // std::cout<<getAngleFromFaceAndVertex(face, edge2(!commun), F)<<"\n";

                    (*harmonic_weights)(v3, 1) = (*harmonic_weights)(v2, 1)
                                                 + ((u(edge1(!commun))-u(edge1(commun)))*cot(getAngleFromFaceAndVertex(face, edge2(last), F)) 
                                                 + (u(edge2(last))-u(edge2(!last)))*cot(getAngleFromFaceAndVertex(face, edge1(!commun), F)))/2.0;
                    n_edges ++;
                    
                }

                else if ((*harmonic_weights)(v3, 1) == -1 && (*harmonic_weights)(v1, 1) != -1){
                    edge1 = E.row(v3);
                    edge2 = E.row(v1);
                    last = 0;
                    if (edge2(last) == edge1(0))  {last = 1; commun = 0;}
                    else if (edge2(last) == edge1(1))  {last = 1; commun = 1;}
                    face = getFaceIndexFromVert(edge1(0), edge1(1), edge2(last), F);
                    // std::cout<<getAngleFromFaceAndVertex(face, edge2(last), F)<<"\n";
                    // std::cout<<getAngleFromFaceAndVertex(face, edge2(!commun), F)<<"\n";

                    (*harmonic_weights)(v3, 1) = (*harmonic_weights)(v1, 1)
                                                 + ((u(edge1(!commun))-u(edge1(commun)))*cot(getAngleFromFaceAndVertex(face, edge2(last), F)) 
                                                 + (u(edge2(last))-u(edge2(!last)))*cot(getAngleFromFaceAndVertex(face, edge1(!commun), F)))/2.0;
                    n_edges ++;
                    
                }

                else if ((*harmonic_weights)(v1, 1) == -1 && (*harmonic_weights)(v3, 1) != -1){
                    edge1 = E.row(v1);
                    edge2 = E.row(v3);
                    last = 0;
                    if (edge2(last) == edge1(0))  {last = 1; commun = 0;}
                    else if (edge2(last) == edge1(1))  {last = 1; commun = 1;}
                    face = getFaceIndexFromVert(edge1(0), edge1(1), edge2(last), F);
                    // std::cout<<getAngleFromFaceAndVertex(face, edge2(last), F)<<"\n";
                    // std::cout<<getAngleFromFaceAndVertex(face, edge2(!commun), F)<<"\n";

                    (*harmonic_weights)(v1, 1) = (*harmonic_weights)(v3, 1)
                                                 + ((u(edge1(!commun))-u(edge1(commun)))*cot(getAngleFromFaceAndVertex(face, edge2(last), F)) 
                                                 + (u(edge2(last))-u(edge2(!last)))*cot(getAngleFromFaceAndVertex(face, edge1(!commun), F)))/2.0;
                    n_edges ++;
                    
                }
            } 
        }

    }

    MatrixXd compute_matrix(MatrixXi &F){

        std::cout<<"Computing Laplacian Matrix...\n";
        MatrixXd L = MatrixXd::Zero(nVertices, nVertices);

        for (int i = 0; i < nVertices; i++){

            float sum = 0;

            MatrixXd adjacency_angles = compute_adjacency(i, F);  //Matrices de taille (n_adj, 3) avec en M(i, )

            for (int j = 0; j < adjacency_angles.rows(); j++){

                int adj_vert = int(adjacency_angles(j, 0));
                if(adjacency_angles(j, 2) == -1.0) {
                    // std::cout<<"Edge at "<<i<<"-"<<adj_vert<<"\n";
                    L(i, adj_vert) = cot(adjacency_angles(j, 1))/2.0;
                    // std::cout<<"Adj : "<<adjacency_angles(j, 1)<<"\n";
                    // std::cout<<"MatrixValue : "<<L(i, adj_vert)<<"\n";
                }
                else {
                    L(i, adj_vert) = (cot(adjacency_angles(j, 1))+cot(adjacency_angles(j, 2)))/2.0;
                    // std::cout<<"Adj : "<<adjacency_angles(j, 1)<<"\n";
                    // std::cout<<"Adj : "<<adjacency_angles(j, 2)<<"\n";
                    // std::cout<<"MatrixValue : "<<L(i, adj_vert)<<"\n";
                }
                sum += L(i, adj_vert);
            }

            L(i, i) = -sum;
            // std::cout<<"MatrixDiagValue : "<<L(i, i)<<"\n";

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
        int res;
		bool v1_in_face, v2_in_face, v3_in_face;

		while (!nf && f < nFaces) {

			v1_in_face = (F.row(f).array() == v1).any();
			v2_in_face = (F.row(f).array() == v2).any();
			v3_in_face = (F.row(f).array() == v3).any();
			if (v1_in_face && v2_in_face && v3_in_face){
				nf = true;
                res = f;
                break;
			}
			f++;
		}

		return res;
	}
    
    int getFaceIndexFromEdge(Vector2i e1, Vector2i e2, MatrixXi &F, int &opp1, int& opp2, int &sommet){
		
		bool nf = false;
		int f = 0;
        int v1 = e1(0);
        int v2 = e1(1);
        int v3;
        if (e2(0) == v1){
            v3 = e2(1);
            opp1 = v2;
            opp2 = v3;
            sommet = v1;
        }
        else if (e2(0) == v2){
            v3 = e2(1);
            opp1 = v1;
            opp2 = v3;
            sommet = v2;

        }
        else if (e2(1) == v1){
            v3 = e2(0);
            opp1 = v2;
            opp2 = v3;
            sommet = v1;
        }
        else{
            v3 = e2(0);
            opp1 = v1;
            opp2 = v3;
            sommet = v2;
        }

        f = getFaceIndexFromVert(v1, v2, v3, F);

        return f;
	}

    double getAngleFromFaceAndVertex(int face, int vertex, MatrixXi &F){
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