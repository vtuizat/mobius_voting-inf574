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
        faces = new MatrixXi(nFaces, 3);

	}

    MidEdgeDS(TriangleMeshDS &mesh)
    {

        nVertices = mesh.sizeOfEdges();
        nFaces = mesh.sizeOfFaces();

        vertices = new MatrixXd(nVertices, 3);
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
        //MatrixXi mesh_edges = mesh.getEdges();
        MatrixXi mesh_faces = mesh.getFaces();
        VectorXi face_v_idx = VectorXi::Zero(nFaces);

        MatrixXi mesh_edges;
        igl::edges(mesh_faces, mesh_edges);

        int vert1_idx, vert2_idx;
        Vector2i face;
        Vector3d mid_vert;

        for (int edge = 0; edge < mesh_edges.rows(); edge++){                
        
            vert1_idx = mesh_edges(edge, 0);
            vert2_idx = mesh_edges(edge, 1);

            face = mesh.getFaceIndexFromEdge(vert1_idx, vert2_idx);

            mid_vert = (mesh_vert.row(vert1_idx)+mesh_vert.row(vert2_idx))/2;

            (*vertices).row(edge) = mid_vert;
            if (face(1) == -1){
                (*faces)(face(0), face_v_idx(face(0))) = edge;
                face_v_idx(face(0))++;
            }
            else {
                (*faces)(face(0), face_v_idx(face(0))) = edge;
                (*faces)(face(1), face_v_idx(face(1))) = edge;
                face_v_idx(face(0))++;
                face_v_idx(face(1))++;
            }

        }

        std::cout<<"All faces completed? "<< (face_v_idx.array() == 3).all()<<"\n";

    }

	int nVertices, nEdges, nFaces;

	MatrixXd *vertices;
	MatrixXi *edges;
	MatrixXi *faces;
};