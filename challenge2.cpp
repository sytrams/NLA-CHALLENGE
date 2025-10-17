#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <cmath>
#include <limits>

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
    cout<< "Frobenius norm of Ag: "<< fb <<endl;

    VectorXd vg= Ag.rowwise().sum();
    //cout<< "vg: "<<"\n" <<vg <<endl;

    //define Dg diagonal matrix with vg as diagonal
    MatrixXd Dg = vg.asDiagonal();
    //cout<< "Dg: "<<"\n" <<Dg <<endl;

    MatrixXd Lg = Dg -Ag;
    VectorXd x = VectorXd::Ones(n);

    // y=Lg*x
    VectorXd y = Lg*x;
    cout<<"Euclidian norm of y: " <<y.norm()<<endl;

    //verify if Lg is SPD
    if (Lg.isApprox(Lg.transpose())) {
        cout << "Lg is symmetric." << endl;
        //control if definite postiive
        LLT<MatrixXd> llt;
        llt.compute(Lg);
        if(llt.info()!= Success) {
            cerr << "Lg isn't positive definite." << endl;
        }else {
            cout << "Lg is positive definite." << endl;
        }
    } else {
        cout << "Lg is not symmetric." << endl;
    }

    SelfAdjointEigenSolver<MatrixXd> eigensolver(Lg); //only computes real eigen values (and Lg is real matrix)
    if (eigensolver.info() != Eigen::Success) abort();
    //cout << "The eigenvalues of Lg:\n" << eigensolver.eigenvalues() << endl;

    VectorXd Lg_eigenvalues = eigensolver.eigenvalues();
    MatrixXd Lg_eigenvectors= eigensolver.eigenvectors();

    //Lg is a real matrix
    cout << "Lg's largest eigenvalue: " << Lg_eigenvalues.maxCoeff() << endl;
    cout << "Lg's smallest eigenvalue: " << Lg_eigenvalues.minCoeff() << endl;

    int smallestPositiveIndex=-1;
    double smallestPositive = numeric_limits<double>::infinity(); //smallestPositive is infinite
    for(int i=0; i<Lg_eigenvalues.size(); i++){
        if (Lg_eigenvalues(i)>0 && Lg_eigenvalues(i)<smallestPositive){
            smallestPositiveIndex=i;
            smallestPositive = Lg_eigenvalues(i);
        }
    }

    if (smallestPositiveIndex == -1) {
        cout << "No positive eigenvalues found." << endl;
    } else {
        cout << "Smallest positive eigenvalue: " << smallestPositive << endl;
        // Corresponding eigenvector: (eigvalue in i position -> eigvector in i column)
        cout << "Corresponding eigenvector:\n" <<  Lg_eigenvectors.col(smallestPositiveIndex) << endl;
    }

    return 0;    
}

