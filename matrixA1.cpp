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

     for (int i = 0; i < height; ++i) {
      for (int j = 0; j < width; ++j) {
        if(smooth(i,j)<0){
          smooth(i,j)=0;}
        if(smooth(i,j)>1){
          smooth(i,j)=1;}
      }
    }

    cout << "Size of A1: " << A1.rows() << "x" << A1.cols() << endl;
    cout << "Number of non-zeros in A1: " << A1.nonZeros() << endl;

    smooth = A1 * w; //smoothing operation

    //add control for values of smooth matrix

     for (int i = 0; i < height; ++i) {
      for (int j = 0; j < width; ++j) {
        if(smooth(i,j)<0){
          smooth(i,j)=0.0;}
        if(smooth(i,j)>1){
          smooth(i,j)=1.0;}
      }
    }
    ....

//now sharpening
     std::vector<T> tripletList2;
    tripletList2.reserve(width*height);
    
      // fill the triplet list with the values for the diagonal matrix
      for( int r =0; r<height; r++){
        for(int c=0; c< width; c++){
          //r*width+c is the index of the pixel in the vectorized form
          tripletList2.push_back(T(r*width+c, r*width+c, 9.0)); //diagonal element
          if(c>0) //left neighbor
            tripletList2.push_back(T(r*width+c, r*width+c-1, -2.0));
          if(c<width-1) //right neighbor
            tripletList2.push_back(T(r*width+c, r*width+c+1,-2.0));
          if(r>0) //top neighbor
            tripletList2.push_back(T(width*(r-1)+c, width*r+c, -2.0));
          if(r<height-1) //bottom neighbor
            tripletList2.push_back(T(width*(r+1)+c, width*r+c, -2.0));
        }
      }
    
    A2.setFromTriplets(tripletList2.begin(), tripletList2.end());

    sharp = A2 * v; //sharping operation


     for (int i = 0; i < height; ++i) {
      for (int j = 0; j < width; ++j) {
        if(sharp(i,j)<0){
          sharp(i,j)=0.0;}
        if(sharp(i,j)>1){
          sharp(i,j)=1.0;}
      }
    }
