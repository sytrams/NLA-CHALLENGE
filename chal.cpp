#include <Eigen/Dense>
#include <iostream>
#include <cstdlib>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include <random>
#include <fstream>
#include <sstream>
#include <string>

// from https://github.com/nothings/stb/tree/master
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include <Eigen/IterativeLinearSolvers>

using namespace Eigen;
using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;

typedef Eigen::Triplet<double> T;

// Function to read a MatrixMarket "vector coordinate real general" file
VectorXd loadCoordinateVector(const string& filename) {
    ifstream file(filename);
    if (!file) {
        throw runtime_error("Could not open file: " + filename);
    }

    string line;

    // Skip comments (lines starting with '%')
    do {
        getline(file, line);
    } while (!line.empty() && line[0] == '%');

    // First non-comment line = vector size
    istringstream iss(line);
    int n;
    iss >> n;

    // Create a dense vector initialized with zeros
    VectorXd vec = VectorXd::Zero(n);

    // Read (index, value) pairs
    int idx;
    double val;
    while (file >> idx >> val) {
        vec(idx - 1) = val;  // MTX uses 1-based indices
    }

    return vec;
}


int main(int argc, char* argv[]) {
  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " <image_path>" << endl;
    return 1;
  }

  const char* input_image_path = argv[1];

  // Load the image using stb_image
  int width, height, channels;
  // for greyscale images force to load only one channel
  unsigned char* image_data = stbi_load(input_image_path, &width, &height, &channels, 1);
  if (!image_data) {
    cerr << "Error: Could not load image " << input_image_path
              << endl;
    return 1;
  }

  cout << "Image loaded: " << height << "x" << width << " with "
            << channels << " channels." << endl;

    // Prepare Eigen matrices for each "modification"
    MatrixXd noise(height, width), mat(height, width), smooth(height, width), sharp(height, width), sobel(height, width), matx(height, width), maty(height, width);
    VectorXd smooth_vec(width*height), sharp_vec(width*height), sobel_vec(width*height);
    //Initialize the noise matrix with values greater than 1 to enter the while loop

    noise.setConstant(2.0);
    //To generate random numbers between -40 and +40
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(-40, +40);

    // Fill the matrices with image data
    for (int i = 0; i < height; ++i) {
      for (int j = 0; j < width; ++j) {
        int index = (i * width + j) * channels;  // 1 channel (Greyscale) 3 channels (RGB)
        mat(i, j) = static_cast<double>(image_data[index]) / 255.0;
        while(noise(i,j)<0 || noise(i,j)>1){
          noise(i,j) = (static_cast<double>(image_data[index])+ dis(gen) )/ 255.0;}
      }
    }

    // Free memory!!!
  stbi_image_free(image_data);

  Matrix<unsigned char, Dynamic, Dynamic, RowMajor> noise_image(height, width);
  // Use Eigen's unaryExpr to map the grayscale values (0.0 to 1.0) to 0 to 255
  noise_image = noise.unaryExpr([](double val) -> unsigned char {
    return static_cast<unsigned char>(val * 255.0);
  });

  // Save the image using stb_image_write
  const string output_image_path1 = "noise_image.png";
  if (stbi_write_png(output_image_path1.c_str(), width, height, 1,
                     noise_image.data(), width) == 0) {
    cerr << "Error: Could not save noise image" << endl;

    return 1;
  }

   cout << "Images saved to " << output_image_path1 << endl;

   VectorXd v(width*height); //image matrix to vector
    for(int i=0; i<height; i++){
      for(int j=0; j<width; j++){
        v(i*width+j)=mat(i,j);
      }
    }
    cout << "m*n: " << width*height << endl;
    cout << "v size: " << v.size() << endl;

    VectorXd w(width*height); //noise image matrix to vector
    for( int i =0; i<height; i++){
      for(int j=0; j<width; j++){
        w(i*width+j)=noise(i,j);
      }
    }
    cout << "w size: " << w.size() << endl;

    cout << "Euclidean norm of v: " << v.norm() << endl;

    SparseMatrix<double> A1(width*height, width*height), A2(width*height, width*height), A3(width*height, width* height);
    vector<T> tripletList;
    tripletList.reserve(width*height);
    
      // fill the triplet list with the values for the diagonal matrix
      for( int r =0; r<height; r++){
        for(int c=0; c< width; c++){
          int index = r*width+c;
          //r*width+c is the index of the pixel in the vectorized form
          tripletList.push_back(T(index, index, 2.0)); //diagonal element
          if(c>0) //left neighbor
            tripletList.push_back(T(index, r*width+c-1, 1.0));
          if(c<width-1) //right neighbor
            tripletList.push_back(T(index, r*width+c+1,1.0));
          if(r>0) //top neighbor
            tripletList.push_back(T(index, width*(r-1)+c, 1.0));
          if(r<height-1) //bottom neighbor
            tripletList.push_back(T(index, width*(r+1)+c, 1.0));
          if(r>0 && c>0) //top-left neighbor
            tripletList.push_back(T(index, width*(r-1)+c-1, 1.0));
          if(r<height-1 && c<width-1) //bottom-right neighbor
            tripletList.push_back(T(index, width*(r+1)+c+1, 1.0));
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
    const string output_image_path2 = "smooth_image.png";
    if (stbi_write_png(output_image_path2.c_str(), width, height, 1,
                     smooth_image.data(), width) == 0) {
      cerr << "Error: Could not save smooth image" << endl;

      return 1;
    }

     cout << "Images saved to " << output_image_path2 << endl;

     vector<T> tripletList2;
    tripletList2.reserve(width*height);
    
      // fill the triplet list with the values for the diagonal matrix
      for( int r =0; r<height; r++){
        for(int c=0; c< width; c++){
          int index =r*width+c;
          //r*width+c is the index of the pixel in the vectorized form
          tripletList2.push_back(T(index, index, 9.0)); //diagonal element
          if(c>0) //left neighbor
            tripletList2.push_back(T(index, r*width+c-1, -2.0));
          if(c<width-1) //right neighbor
            tripletList2.push_back(T(index, r*width+c+1,-2.0));
          if(r>0) //top neighbor
            tripletList2.push_back(T(index, width*(r-1)+c, -2.0));
          if(r<height-1) //bottom neighbor
            tripletList2.push_back(T(index, width*(r+1)+c, -2.0));
        }
      }
    
    A2.setFromTriplets(tripletList2.begin(), tripletList2.end());

    cout << "Number of non-zeros in A2: " << A2.nonZeros() << endl;

    if (A2.isApprox(A2.transpose())) {
    cout << "A2 is symmetric." << endl;
    //control if definite postiive
    SimplicialLDLT<SparseMatrix<double>> solver;
    solver.compute(A2);
    if(solver.info() != Success) {
        cerr << "Decomposition failed, matrix might not be positive definite." << endl;
    } else {
        cout << "Matrix is positive definite." << endl;
    }

} else {
    cout << "A2 is not symmetric." << endl;
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
    const string output_image_path3 = "sharp_image.png";
    if (stbi_write_png(output_image_path3.c_str(), width, height, 1,
                     sharp_image.data(), width) == 0) {
      cerr << "Error: Could not save sharp image" << endl;

      return 1;
    }

     cout << "Images saved to " << output_image_path3 << endl;

    //Export A2 and w in the .mtx format
    

   // Export vector w in .mtx format
    int n = w.size();
    // saveMarketVector(w, "./w.mtx");
    FILE* out = fopen("w.mtx","w");
    fprintf(out,"%%%%MatrixMarket vector coordinate real general\n");
    fprintf(out,"%d\n", n);
    for (int i=0; i<n; i++) {
        fprintf(out,"%d %f\n", i ,w(i));
    }
    fclose(out);

    if (!saveMarket(A2, "A2.mtx")) {
        cerr << "Error: Could not save matrix A2 to A2.mtx" << endl;
        return 1;
    }

       // Load x.mtx into a dense VectorXd
        VectorXd x = loadCoordinateVector("x.mtx");

        cout << "Loaded vector x of size: " << x.size() << endl;
        // x is now available in memory as a VectorXd

    for(int i=0; i<height; i++){
      for(int j=0; j<width; j++){
        matx(i,j)=x(i*width+j);
        if(matx(i,j)<0) {matx(i,j)=0.0;}
        if(matx(i,j)>1) {matx(i,j)=1.0;}
      }
    }
    //cout << matx <<endl;
    
    //VectorXd w_hat = A2*x;
     // cout<< w_hat <<endl;

      //cout << w <<endl;

      //cout<< "Norm: "<< (w_hat-w).norm()<< endl;

     Matrix<unsigned char, Dynamic, Dynamic, RowMajor> matx_image(height, width);
  // Use Eigen's unaryExpr to map the grayscale values (0.0 to 1.0) to 0 to 255
    matx_image = matx.unaryExpr([](double val) -> unsigned char {
    return static_cast<unsigned char>(val * 255.0);
    });

    //cout << "matx: "<< endl << matx << endl;
  // Save the image using stb_image_write
    const string output_image_path5 = "matx_image.png";
    if (stbi_write_png(output_image_path5.c_str(), width, height, 1,
                     matx_image.data(), width) == 0) {
      cerr << "Error: Could not save matx image" << endl;

      return 1;
    }

     cout << "Images saved to " << output_image_path5 << endl;


      vector<T> tripletList3;
    tripletList3.reserve(width*height);
    
      // fill the triplet list with the values for the diagonal matrix
      for( int r =0; r<height; r++){
        for(int c=0; c< width; c++){
          int index= r*width+c;
          //r*width+c is the index of the pixel in the vectorized form
          if(r>0) //top neighbor
            tripletList3.push_back(T(index, width*(r-1)+c, -2.0));
          if(r<height-1) //bottom neighbor
            tripletList3.push_back(T(index, width*(r+1)+c, 2.0));
          if(r>0 && c>0) //top-left neighbor
            tripletList3.push_back(T(index, width*(r-1)+c-1, -1.0));
          if(r<height-1 && c<width-1) //bottom-right neighbor
            tripletList3.push_back(T(index, width*(r+1)+c+1, 1.0));
          if(r>0 && c<width-1) //top-right neighbor
            tripletList3.push_back(T(index, width*(r-1)+c+1, -1.0));
          if(r<height-1 && c>0) //bottom-left neighbor
            tripletList3.push_back(T(index, width*(r+1)+c-1, 1.0));
        }
      }
    
    A3.setFromTriplets(tripletList3.begin(), tripletList3.end());

    if (A3.isApprox(A3.transpose())) {
    cout << "A3 is symmetric." << endl;
} else {
    cout << "A3 is not symmetric." << endl;
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
    const string output_image_path4 = "sobel_image.png";
    if (stbi_write_png(output_image_path4.c_str(), width, height, 1,
                     sobel_image.data(), width) == 0) {
      cerr << "Error: Could not save sobel image" << endl;

      return 1;
    }

     cout << "Images saved to " << output_image_path4 << endl;

    VectorXd y(A3.rows()), yj(A3.rows());
   int result;

   //Identity matrix
  SparseMatrix<double> I(A3.rows(), A3.cols());
  I.setIdentity();

  //Construct BiCG solver

  /*BiCGSTAB<SparseMatrix<double>, IncompleteLUT<double>> BiCG;
  BiCG.setTolerance(1.e-8);
  BiCG.setMaxIterations(5000);
  BiCG.compute(3*I+A3);
  
  y = BiCG.solve(w);

  cout<<"Eigen native BiCG with ILU"<<endl;
  cout<< "#iterations:  " <<BiCG.iterations()<<endl;
  cout<< "estimated error:  "<< BiCG.error()<<endl;
  cout<<"effective error:  "<< (w-y).norm() <<endl;
*/
  BiCGSTAB<SparseMatrix<double>,  DiagonalPreconditioner<double> > BiCGj;
  BiCGj.setTolerance(1.e-8);
  BiCGj.setMaxIterations(5000);
  BiCGj.compute(3*I+A3);
  
  yj = BiCGj.solve(w);

  cout<<"Eigen native BiCG with jacobi"<<endl;
  cout<< "#iterations:  " <<BiCGj.iterations()<<endl;
  cout<< "estimated error:  "<< BiCGj.error()<<endl;
  cout<<"effective error:  "<< (w-yj).norm() <<endl;
  
  for(int i=0; i<height; i++){
      for(int j=0; j<width; j++){
        maty(i,j)=yj(i*width+j);
        if(maty(i,j)<0) {maty(i,j)=0.0;}
        if(maty(i,j)>1) {maty(i,j)=1.0;}
      }
    }
  
  Matrix<unsigned char, Dynamic, Dynamic, RowMajor> maty_image(height, width);
  // Use Eigen's unaryExpr to map the grayscale values (0.0 to 1.0) to 0 to 255
    maty_image = maty.unaryExpr([](double val) -> unsigned char {
    return static_cast<unsigned char>(val * 255.0);
    });

  // Save the image using stb_image_write
    const string output_image_path6 = "maty_image.png";
    if (stbi_write_png(output_image_path6.c_str(), width, height, 1,
                     maty_image.data(), width) == 0) {
      cerr << "Error: Could not save maty image" << endl;

      return 1;
    }

     cout << "Images saved to " << output_image_path6 << endl;

    
  return 0;

}
