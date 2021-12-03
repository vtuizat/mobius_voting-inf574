#include <igl/opengl/glfw/Viewer.h>
#include <math.h>
#include <ostream>


using namespace Eigen;

/**
 * A class for representing linear transformations on 3D points (using homogeneous coordinates)
 * */
class MeshTransformation
{
public:
/*
Initialize the identity transformation
**/
  MeshTransformation()
  {
    MatrixXd m(4, 4);
    m(0, 0) = 1.0; m(1, 0) = 0.0; m(2, 0) = 0.0; m(3, 0) = 0.0;
    m(0, 1) = 0.0; m(1, 1) = 1.0; m(2, 1) = 0.0; m(3, 1) = 0.0;
    m(0, 2) = 0.0; m(1, 2) = 0.0; m(2, 2) = 1.0; m(3, 2) = 0.0;
    m(0, 3) = 0.0; m(1, 3) = 0.0; m(2, 3) = 0.0; m(3, 3) = 1.0;

    M = m;
  }


/*
Initialize a scaling transformation
**/
  MeshTransformation(double s1, double s2, double s3)
  {
    // TO BE COMPLETED
    MatrixXd m(4, 4);
    m(0, 0) = s1 ; m(1, 0) = 0.0; m(2, 0) = 0.0; m(3, 0) = 0.0;
    m(0, 1) = 0.0; m(1, 1) = s2 ; m(2, 1) = 0.0; m(3, 1) = 0.0;
    m(0, 2) = 0.0; m(1, 2) = 0.0; m(2, 2) = s3 ; m(3, 2) = 0.0;
    m(0, 3) = 0.0; m(1, 3) = 0.0; m(2, 3) = 0.0; m(3, 3) = 1.0;

    M = m;
  }

/*
Initialize a rotation transformation around a given axis (X, Y or Z) <br><br>

 @param  direction  a value 0, 1 or 2 indicating the direction (X, Y or Z respectively)
**/
  MeshTransformation(double theta, int direction)
  {
// TO BE COMPLETED
    MatrixXd m(4, 4);

    double cos_theta = std::cos(theta);
    double sin_theta = std::sin(theta);

    if(direction==0){
        m(0, 0) = 1.0; m(1, 0) = 0.0; m(2, 0) = 0.0; m(3, 0) = 0.0;
        m(0, 1) = 0.0; m(1, 1) = cos_theta; m(2, 1) = -1 * sin_theta; m(3, 1) = 0.0;
        m(0, 2) = 0.0; m(1, 2) = sin_theta; m(2, 2) = cos_theta; m(3, 2) = 0.0;
        m(0, 3) = 0.0; m(1, 3) = 0.0; m(2, 3) = 0.0; m(3, 3) = 1.0;
    }

    if(direction==1){
        m(0, 0) = cos_theta; m(1, 0) = 0.0; m(2, 0) = sin_theta; m(3, 0) = 0.0;
        m(0, 1) = 0.0; m(1, 1) = 1.0; m(2, 1) = 0.0; m(3, 1) = 0.0;
        m(0, 2) = -1 * sin_theta; m(1, 2) = 0.0; m(2, 2) = cos_theta; m(3, 2) = 0.0;
        m(0, 3) = 0.0; m(1, 3) = 0.0; m(2, 3) = 0.0; m(3, 3) = 1.0;
    }

    if(direction==2){
        m(0, 0) = cos_theta; m(1, 0) = -1 * sin_theta; m(2, 0) = 0.0; m(3, 0) = 0.0;
        m(0, 1) = sin_theta; m(1, 1) = cos_theta; m(2, 1) = 0.0; m(3, 1) = 0.0;
        m(0, 2) = 0.0; m(1, 2) = 0.0; m(2, 2) = 1.0; m(3, 2) = 0.0;
        m(0, 3) = 0.0; m(1, 3) = 0.0; m(2, 3) = 0.0; m(3, 3) = 1.0;
    }

    M = m;
  }

/*
Initialize a translation
**/
  MeshTransformation(RowVector3d t)
  {
// TO BE COMPLETED
    MatrixXd m(4, 4);
    m(0, 0) = 1.0; m(1, 0) = 0.0; m(2, 0) = 0.0; m(3, 0) = 0.0;
    m(0, 1) = 0.0; m(1, 1) = 1.0; m(2, 1) = 0.0; m(3, 1) = 0.0;
    m(0, 2) = 0.0; m(1, 2) = 0.0; m(2, 2) = 1.0; m(3, 2) = 0.0;
    m(0, 3) = t(0); m(1, 3) = t(1); m(2, 3) = t(2); m(3, 3) = 1.0;

    M = m;
  }

/*
Matrix accessor

@return  the matrix transformation
**/
  MatrixXd get_matrix() {
    return M;
  }

/*
Initialize a transformation given an input matrix 'm'
**/
  void set_matrix(MatrixXd m)
  {
    M = m;
  }

/*
Apply the transformation to all vertices stored in a matrix 'V' <br>

@param V  vector storing the input points
**/
  void transform(MatrixXd &V) {
    // TO BE COMPLETED
    for (unsigned i = 0; i < V.rows(); i++){
      V.row(i)=transform(V.row(i));
    }
  }

  	/**
	 * Apply the transformation to a 3d (row) vector 'v' <br>
   * 
   * Remark: use homogeneous coordinates
   * 
   * @return  the vector after transformation
	 */
	RowVector3d transform(RowVector3d v) {
    // TO BE COMPLETED
    RowVector3d result(0., 0., 0.);
    RowVector4d hResult(v(0), v(1), v(2), 1.);
    hResult = M * hResult.transpose();
    result << hResult(0), hResult(1), hResult(2);
		return result;
	}

	/**
	 * Compose the current transformation with a transfomation 't': return a new transformation
	 */
	MeshTransformation compose(MeshTransformation t) {
    MeshTransformation res(1.0, 1.0, 1.0);
    // TO BE COMPLETED
    res.set_matrix(M * t.get_matrix());
    return res;
	}

	/**
	 * Print the matrix transformation
	 */
  friend std::ostream& operator<<(std::ostream &os, MeshTransformation& t) {
    return os << "matrix:\n" << t.get_matrix() << std::endl;
  }

private:
  MatrixXd M; // a 4x4 matrix representing a linear transformation
};
