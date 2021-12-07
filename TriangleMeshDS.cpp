#include <igl/opengl/glfw/Viewer.h>
#include <igl/edges.h>
#include <igl/per_vertex_normals.h>
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

		N = new MatrixXd(nVertices, 3);
		C = new MatrixXd(nVertices, 3);

        if (compute_normals) {
            compute_per_vertex_normals();
        }

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

	MatrixXd getVert(){
		return *V;
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
		float cosa = v1.dot(v2);
		if(cosa >= 1.f)
			return 0.f;
		else if(cosa <= -1.f)
			return M_PI;
		else
			return std::acos(cosa);
	}

	int nVertices, nFaces;

	MatrixXd *V;
	MatrixXi *F;
    MatrixXd *N;
    MatrixXd *C;

};