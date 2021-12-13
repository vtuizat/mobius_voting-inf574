#include <igl/opengl/glfw/Viewer.h>
#include <igl/edges.h>
#include <igl/per_vertex_normals.h>
#include <igl/exact_geodesic.h>
#include <igl/edges.h>
#include <fstream>
#include <iostream>
#include <ostream>
#include <vector>
#include <cassert>


using namespace Eigen;
using namespace std;

class TriangleMeshDS
{

public:

	/** 
	 * Allocate the memory space for storing the data structure.
	 * This version of the constructor stores face/halfedge incidence relations
	 **/
	TriangleMeshDS(int n, int f, bool compute_normals = false)
	{
		nVertices = n;
		nFaces = f;
        
		V = new MatrixXd(nVertices, 3);
		F = new MatrixXi(nFaces, 3);
		N = new MatrixXd(nVertices, 3);
		C = new MatrixXd(nVertices, 3);

        if (compute_normals) {
            compute_per_vertex_normals();
        }

	}

    TriangleMeshDS(MatrixXd &V_in, MatrixXi &F_in, bool compute_normals = false)
    {
        nVertices = V_in.rows();
        nFaces = F_in.rows();

        V = &V_in;
    	F = &F_in;

		igl::edges(*F, *E);
		std::cout<<(*E).rows()<<" "<<(*E).cols()<<"\n";

		nEdges = (*E).rows();

		std::cout<<" Number of edges : "<<nEdges<<"\n";

		N = new MatrixXd(nVertices, 3);
		C = new MatrixXd(nVertices, 3);

        if (compute_normals) {
            compute_per_vertex_normals();
        }

		nVertices = (*V).rows();

    }

	MatrixXd compute_angles(){

		std::cout<<"Computing angles...\n";


		MatrixXd angles_per_faces(nFaces, 3);
		
		float angle;

		for(int i = 0; i < nFaces; i++){
			for(int idx = 0; idx < 3; idx++){
				
				int opp1 = (idx+1)%3;
				int opp2 = (idx+2)%3;
				
				int vert_idx = (*F)(i,idx);
				int vert_opp1 = (*F)(i, opp1);
				int vert_opp2 = (*F)(i, opp2);

				angle = angle_between((*V).row(vert_opp1)-(*V).row(vert_idx), (*V).row(vert_opp2)-(*V).row(vert_idx));

				angles_per_faces(i, idx) = angle;
			}
		}

		return angles_per_faces;

	}

	int get_cut_face(){

		std::cout<<"Compute Face to cut off...\n";

		double mean_dist, min_mean;

		VectorXd dist_f;
		VectorXi all_verts;
		all_verts.setLinSpaced((*V).rows(),0,(*V).rows()-1);
		dist_f = geodesic_dist_face(0, all_verts);
		min_mean = dist_f.mean();

		int idx = 0;
		for (int v = 1; v < nVertices; v++){
			dist_f = geodesic_dist_face(v, all_verts);
			mean_dist = dist_f.mean();
			if (mean_dist < min_mean) {
				min_mean = mean_dist;
				idx = v;
			}
		}

		int face_idx = 0;
		while ((*F)(face_idx, 0)!=idx && (*F)(face_idx, 1)!=idx && (*F)(face_idx, 1)!=idx){
			face_idx++;
		}

		return face_idx;
	}

	int get_random_cut_off_face(){
		int face = std::rand()%nFaces;
		std::cout<<"Selected face to cut off : "<<face<<"\n";
		return face;
	}

	VectorXd geodesic_dist_face(int v, VectorXi targets){

		VectorXi VS, FS, VT, FT;
		VS.resize(1);
		VS << v;
		VT = targets;
		VectorXd dist;
		igl::exact_geodesic(*V, *F, VS, FS, VT, FT, dist);

		return dist;

	}

	Vector2i getFaceIndexFromEdge(int v1, int v2){
		
		Vector2i idx = {-1, -1};
		int nf = 0;
		int f = 0;
		bool v1_in_face, v2_in_face;

		while (nf < 2 && f < nFaces) {

			v1_in_face = ((*F).row(f).array() == v1).any();
			v2_in_face = ((*F).row(f).array() == v2).any();
			if (v1_in_face && v2_in_face){
				idx(nf) = f;
				nf++;
			}
			f++;
		}

		return idx;
	}

	// int getFacefromEdge(int edge){
	// }

	MatrixXd getVert(){
		return *V;
	}	
	
	MatrixXi getEdges(){
		return *E;
	}

	MatrixXi getFaces(){
		return *F;
	}

	/** 
	 * Return the number of vertices
	 **/
	int sizeOfVertices()
	{
		return nVertices;
	}

	/** 
	 * Return the number of vertices
	 **/
	int sizeOfEdges()
	{
		return nEdges;
	}


	/** 
	 * Return the number of faces
	 **/
	int sizeOfFaces()
	{
		return nFaces;
	}

	/** 
	 * Print the array T[] storing all references
	 **/
	void print()
	{
        
            std::cout<<"Number of vertices : "<<nVertices<<"\n";
            std::cout<<"Number of faces : "<<nFaces<<"\n\n";

	}

private:

    void compute_per_vertex_normals(){

		std::cout<<"Computing normals per vertex...\n";

        igl::per_vertex_normals(*V, *F, *N);
	
    }

	float angle_between(Vector3d v1, Vector3d v2)
	{
		float dot = v1(0)*v2(0) + v1(1)*v2(1) + v1(2)*v2(2);    // #between [x1, y1, z1] and [x2, y2, z2]
		float lenSq1 = v1(0)*v1(0) + v1(1)*v1(1) + v1(2)*v1(2);
		float lenSq2 = v2(0)*v2(0) + v2(1)*v2(1) + v2(2)*v2(2);
		float angle = acos(dot/std::sqrt(lenSq1 * lenSq2));

		return angle;
	}

	int nVertices, nEdges, nFaces;

	MatrixXd *V;
	MatrixXi *E;
	MatrixXi *F;
    MatrixXd *N;
    MatrixXd *C;

};