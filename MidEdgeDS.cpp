#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <ostream>

#ifndef TRIANGLEMESHDS_HEADER
   #define TRIANGLEMESHDS_HEADER
   #include "TriangleMeshDS.cpp"
#endif

using namespace Eigen;
using namespace std;

class MidEdgeDS
{

public:

	/** 
	 * Allocate the memory space for storing the data structure.
	 * This version of the constructor stores face/halfedge incidence relations
	 **/
	MidEdgeDS(int n, int e, int f)
	{
		nVertices = n;
		nEdges = e;
		nFaces = f;
        
        vertices = new MatrixXd(nVertices, 3);
        //edges = new MatrixXi(nEdges, 2);
        faces = new MatrixXi(nFaces, 3);

	}

    MidEdgeDS(TriangleMeshDS &mesh)
    {

        nVertices = (mesh.sizeOfVertices()-1)*3;
        nEdges = mesh.sizeOfFaces()*3;
        nFaces = mesh.sizeOfFaces();

        vertices = new MatrixXd(nVertices, 3);
        // edges = new MatrixXi(nEdges, 2);
        faces = new MatrixXi(nFaces, 3);

        computeMidEdgeDS(mesh);
    }



	//--- methods for accessing data and navigating between mesh elements ---

	MatrixXd getVert(){
		return *vertices;
	}

	MatrixXi getFaces(){
		return *faces;
	}


	/** 
	 * Return the number of vertices
	 **/
	int sizeOfVertices()
	{
		return nVertices;
	}

	/** 
	 * Return the number of half-edges
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
            std::cout<<"Number of edges : "<<nEdges<<"\n";
            std::cout<<"Number of faces : "<<nFaces<<"\n";

	}

private:

    void computeMidEdgeDS(TriangleMeshDS &mesh){

        std::cout<<"Computing Mid Edge Mesh...\n";

        MatrixXd mesh_vert = mesh.getVert();
        MatrixXi mesh_faces = mesh.getFaces();

        MatrixXi check_edge = MatrixXi::Constant(nEdges, nEdges, -1);

        int count_vertices = 0;

        for (int face = 0; face < nFaces; face++){

            Vector3i vert_idxs = mesh_faces.row(face);

            Vector3i new_vert_idxs;
            // Vector2i new_edge1, new_edge2, new_edge3;
            Vector3i new_face;

            for(int edge_in_face = 0; edge_in_face < 3; edge_in_face++){

                int vert1_idx = vert_idxs(edge_in_face);
                int vert2_idx = vert_idxs((edge_in_face+1)%3);

                int new_vert_idx;

                if (check_edge(vert1_idx, vert2_idx) == -1){
                    
                    Vector3d mid_vert;

                    Vector3d vert1 = mesh_vert.row(vert1_idx);
                    Vector3d vert2 = mesh_vert.row(vert2_idx);

                    mid_vert = (vert1+vert2)/2;

                    (*vertices).row(count_vertices) = mid_vert;

                    check_edge(vert1_idx, vert2_idx) = count_vertices;
                    check_edge(vert2_idx, vert1_idx) = count_vertices;
                    count_vertices++;
                }
                
                new_vert_idx = check_edge(vert1_idx, vert2_idx);

                new_vert_idxs(edge_in_face) = new_vert_idx;

            }

            // new_edge1(0) = new_vert_idxs(0);
            // new_edge1(1) = new_vert_idxs(1);

            
            // new_edge2(0) = new_vert_idxs(1);
            // new_edge2(1) = new_vert_idxs(2);

            // new_edge3(0) = new_vert_idxs(2);
            // new_edge3(1) = new_vert_idxs(0);

            // (*edges).row(face*3) = new_edge1;
            // (*edges).row(face*3+1) = new_edge2;
            // (*edges).row(face*3+2) = new_edge3;

            new_face = {new_vert_idxs(0), new_vert_idxs(1), new_vert_idxs(2)};
            (*faces).row(face) = new_face;

        }

    }

	int nVertices, nEdges, nFaces;

	MatrixXd *vertices;
	//MatrixXi *edges;
	MatrixXi *faces;
};