#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <cmath>

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{
    //create adjacency matrix
    int n = 9;	
    SparseMatrix<double> mat(9, 9);                           // define matrix
    mat.coeffRef(0,1)= 1.0;
    mat.coeffRef(0,3)=1.0;
    mat.coeffRef(1, 2)=1.0;
    mat.coeffRef(2, 3)=1.0;
    mat.coeffRef(2, 4)=1.0;
    mat.coeffRef(4, 5)=1.0;
    mat.coeffRef(4, 7)=1.0;
    mat.coeffRef(4,8)= 1.0;
    mat.coeffRef(5, 6)=1.0;
    mat.coeffRef(6, 7)=1.0;
    mat.coeffRef(6, 8)=1.0;
    mat.coeffRef(7, 8)=1.0;

    MatrixXd A, Ag;
    A = MatrixXd(mat);
    Ag= A + A.transpose(); //define Adj matrix
    //cout << Ag << endl;

    //Frobenius norm, Ag is symmetric
    double fb = sqrt((Ag*Ag).trace());
    cout<< "Frobenius norm: "<< fb <<endl;


    return 0;    
}
