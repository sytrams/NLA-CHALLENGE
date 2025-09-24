#include <iostream>
#include <cstdlib>
#include <Eigen/Dense>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

using namespace Eigen;

using Eigen::MatrixXd;
using Eigen::VectorXd;

int main(int argc, char* argv[]) 
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <image_path>" << std::endl;
    return 1;
  }

  const char* input_image_path = argv[1];

  // Load the image using stb_image
  int width, height, channels;
  unsigned char* image_data = stbi_load(input_image_path, &width, &height, &channels, 1);  // Force load as RGB

  if (!image_data) {
    std::cerr << "Error: Could not load image " << input_image_path << std::endl;
    return 1;
  }

  //Report size of the matrix
  std::cout << "Image loaded: " << width << "x" << height << " with " << channels << " channels." << std::endl;

  //Prepare Eigen matrices
  MatrixXd noise(height, width), mat(height, width);

  VectorXd w(height*width), v(height*width);

  //Fill matrices
  for (int i=0; i<height; i++)
  {
    for (int j=0; j<width; j++)
    {
        int index = (i*width+j) * channels; //1 channel greyscale, 3 channels RGB
        int jump=0;
        mat(i,j)=static_cast<double>(image_data[index]);  //copy image in matrix
        noise(i,j)= mat(i,j) + (rand() % 101) - 40;  // Add causal noise [-40,40] to the matrix
        v(index) = mat(i,j);
        w(index) = noise(i,j);
        if (noise(i,j)<0.0)
        {
            noise(i,j)=0.0;
            w(index)=0.0;
        }
        if (noise(i,j)>255.0)
        {
            noise(i,j)=255.0;
            w(index)=255.0;
        }
    }
  }

  //Verify the size of the vectors
  std::cout << "Vector v has size: " << v.size() << std::endl;
  std::cout << "Vector w has size: " << w.size() << std::endl;

  //Euclidean norm of v
  std::cout << "Euclidean norm of v: " << v.norm() << std::endl;

  // Free memory!!!
  stbi_image_free(image_data);

  Matrix<unsigned char, Dynamic, Dynamic, RowMajor> noise_image(height,width);

  //use Eigen's unaryExpr to map the greyscale values (0.0 to 1.0) to 0 to 255
  noise_image = noise.unaryExpr([](double val) -> unsigned char
  {
    return static_cast<unsigned char>(val);
  });
  
  //Save the image using stb_image_write
  const std::string output_image_path1 = "noise_image.png";
  if (stbi_write_png(output_image_path1.c_str(),width,height,1,noise_image.data(),width)==0)
  {
    std::cout << "Error: Could not save greyscale image" << std::endl;

    return 1;
  }
  
  return 0;
}
