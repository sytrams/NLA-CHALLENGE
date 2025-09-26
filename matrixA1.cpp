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

 // Prepare Eigen matrices for each "modification"
    MatrixXd noise(height, width), mat(height, width), smooth(height, width), sharp(height, width), sobel(height, width);
    VectorXd smooth_vec(width*height), sharp_vec(width*height), sobel_vec(width*height);
    //Initialize the noise matrix with values greater than 1 to enter the while loop

    SparseMatrix<double> A1(width*height, width*height), A2(width*height, width*height), A3(width*height, width* height);
    
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

    smooth_vec = A1 * w; //smoothing operation

    for(int i=0; i<height; i++){
      for(int j=0; j<width; j++){
        smooth(i,j)=smooth_vec(i*width+j);
        if(smooth(i,j)<0) {smooth(i,j)=0.0;}
        if(smooth(i,j)>1) {smooth(i,j)=1.0;}
      }
    }
    

   Matrix<unsigned char, Dynamic, Dynamic, RowMajor> smooth_image(height, width);
  // Use Eigen's unaryExpr to map the grayscale values (0.0 to 1.0) to 0 to 255
    smooth_image = smooth.unaryExpr([](double val) -> unsigned char {
    return static_cast<unsigned char>(val * 255.0);
    });

  // Save the image using stb_image_write
    const std::string output_image_path2 = "smooth_image.png";
    if (stbi_write_png(output_image_path2.c_str(), width, height, 1,
                     smooth_image.data(), width) == 0) {
      std::cerr << "Error: Could not save smooth image" << std::endl;

      return 1;
    }

     std::cout << "Images saved to " << output_image_path2 << std::endl;

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
    cout << "Number of non-zeros in A2: " << A2.nonZeros() << endl;

    if (A2.isApprox(A2.transpose())) {
     std::cout << "A2 is symmetric." << std::endl;
    } else {
    std::cout << "A2 is not symmetric." << std::endl;
    }

    sharp_vec = A2 * v; //sharping operation

    for(int i=0; i<height; i++){
      for(int j=0; j<width; j++){
        sharp(i,j)=sharp_vec(i*width+j);
        if(sharp(i,j)<0) {sharp(i,j)=0.0;}
        if(sharp(i,j)>1) {sharp(i,j)=1.0;}
      }
    }
    

   Matrix<unsigned char, Dynamic, Dynamic, RowMajor> sharp_image(height, width);
  // Use Eigen's unaryExpr to map the grayscale values (0.0 to 1.0) to 0 to 255
    sharp_image = sharp.unaryExpr([](double val) -> unsigned char {
    return static_cast<unsigned char>(val * 255.0);
    });

  // Save the image using stb_image_write
    const std::string output_image_path3 = "sharp_image.png";
    if (stbi_write_png(output_image_path3.c_str(), width, height, 1,
                     sharp_image.data(), width) == 0) {
      std::cerr << "Error: Could not save sharp image" << std::endl;

      return 1;
    }

     std::cout << "Images saved to " << output_image_path3 << std::endl;

      std::vector<T> tripletList3;
    tripletList3.reserve(width*height);
    
      // fill the triplet list with the values for the diagonal matrix
      for( int r =0; r<height; r++){
        for(int c=0; c< width; c++){
          //r*width+c is the index of the pixel in the vectorized form
          if(r>0) //top neighbor
            tripletList3.push_back(T(width*(r-1)+c, width*r+c, -2.0));
          if(r<height-1) //bottom neighbor
            tripletList3.push_back(T(width*(r+1)+c, width*r+c, 2.0));
          if(r>0 && c>0) //top-left neighbor
            tripletList3.push_back(T(width*(r-1)+c-1, width*r+c, -1.0));
          if(r<height-1 && c<width-1) //bottom-rigtht neighbor
            tripletList3.push_back(T(width*(r+1)+c+1, width*r+c, 1.0));
          if(r>0 && c<width-1) //top-right neighbor
            tripletList3.push_back(T(width*(r-1)+c+1, width*r+c, -1.0));
          if(r<height-1 && c>0) //bottom-left neighbor
            tripletList3.push_back(T(width*(r+1)+c-1, width*r+c, 1.0));
        }
      }
    
    A3.setFromTriplets(tripletList3.begin(), tripletList3.end());

    if (A3.isApprox(A3.transpose())) {
     std::cout << "A3 is symmetric." << std::endl;
    } else {
     std::cout << "A3 is not symmetric." << std::endl;
    }

    sobel_vec = A3 * v; //sharping operation

    for(int i=0; i<height; i++){
      for(int j=0; j<width; j++){
        sobel(i,j)=sobel_vec(i*width+j);
        if(sobel(i,j)<0) {sobel(i,j)=0.0;}
        if(sobel(i,j)>1) {sobel(i,j)=1.0;}
      }
    }
    

   Matrix<unsigned char, Dynamic, Dynamic, RowMajor> sobel_image(height, width);
  // Use Eigen's unaryExpr to map the grayscale values (0.0 to 1.0) to 0 to 255
    sobel_image = sobel.unaryExpr([](double val) -> unsigned char {
    return static_cast<unsigned char>(val * 255.0);
    });

  // Save the image using stb_image_write
    const std::string output_image_path4 = "sobel_image.png";
    if (stbi_write_png(output_image_path4.c_str(), width, height, 1,
                     sobel_image.data(), width) == 0) {
      std::cerr << "Error: Could not save sobel image" << std::endl;

      return 1;
    }

     std::cout << "Images saved to " << output_image_path4 << std::endl;



  return 0;
}
