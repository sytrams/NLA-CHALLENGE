#include <Eigen/Dense>
#include <iostream>
#include <cstdlib>

#include <unsupported/Eigen/SparseExtra>
#include <random>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

using namespace Eigen;
using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;

typedef Eigen::Triplet<double> T;

//da mettere dentro l'int main
int main(...):{

SparseMatrix<double> A1(width*height, width*height);
    
    std::vector<T> tripletList;
    tripletList.reserve(width*height);
    
      // fill the triplet list with the values for the diagonal matrix
      for( int r =0; r<height; r++){
        for(int c=0; c< width; c++){
          //r*width+c is the index of the pixel in the vectorized form
          tripletList.push_back(T(r*width+c, r*width+c, 2.0)); //diagonal element
          if(c>0) //left neighbor
            tripletList.push_back(T(r*width+c, r*width+c-1, 1.0));
          if(c<width-1) //right neighbor
            tripletList.push_back(T(r*width+c, r*width+c+1,1.0));
          if(r>0) //top neighbor
            tripletList.push_back(T(width*(r-1)+c, width*r+c, 1.0));
          if(r<height-1) //bottom neighbor
            tripletList.push_back(T(width*(r+1)+c, width*r+c, 1.0));
          if(r>0 && c>0) //top-left neighbor
            tripletList.push_back(T(width*(r-1)+c, width*r+c-1, 1.0));
          if(r<height-1 && c<width-1) //bottom-rigtht neighbor
            tripletList.push_back(T(width*(r+1)+c, width*r+c+1, 1.0));
        }
      }
    
    A1.setFromTriplets(tripletList.begin(), tripletList.end());
    A1=A1/8.0;

    cout << "Size of A1: " << A1.rows() << "x" << A1.cols() << endl;
    cout << "Number of non-zeros in A1: " << A1.nonZeros() << endl;

    smooth = A1 * w; //smoothing operation
