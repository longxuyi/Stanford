#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>
#include <iostream>
#include <jpeglib.h>
#include <sstream>
#include <stdexcept>
#include <string>

#include "image.hpp"
#include "hw6.hpp"

Image::Image(std::string filename) {
    
    this->input_file = filename;
    ReadGrayscaleJPEG(this->input_file, this->img);
    this->input.resize(boost::extents[this->img.shape()[0]][this->img.shape()[1]]);
    this->output.resize(boost::extents[this->img.shape()[0]][this->img.shape()[1]]);
    this->input = this->img;
    this->output = this->input;
}

void Image::Save(std::string filename) {
    if (filename.empty()) {
        filename = this->input_file;
    }
    WriteGrayscaleJPEG(filename, this->output);
}

void Image::Convolution(boost::multi_array<unsigned char,2>& input, 
                        boost::multi_array<unsigned char,2>& output,
                        boost::multi_array<float,2>& kernel) {

    int nRow = (int) input.shape()[0];
    int nCol = (int) input.shape()[1];
    int nKer = (int) kernel.shape()[0];

    if(nRow!=(int)output.shape()[0] || nCol!=(int)output.shape()[1]) {
        std::cerr << "error: input and output are not the same size" << std::endl;
        exit(1);
    }

    if(nKer%2 == 0 || nKer <3 || kernel.shape()[0] != kernel.shape()[1]) {
        std::cerr << "error: kernel size must be square, odd size, and at least size 3" << std::endl;
        exit(1);        
    }
    
    // For each element in the input array
    //determine value for each element in the output array
    
    for (int i = 0; i < nRow; i++) {
        for (int  j = 0; j < nCol; j++) {
            //convolution calculation
            float temp = 0.f;
            int kr=0, kc =0;

            for (int m = i - (nKer-1)/2; m <= i + (nKer-1)/2; m++) {
                int r = m;
                if (r < 0) r = 0;
                if (r > nRow-1) r = nRow-1;

                for (int n = j - (nKer-1)/2; n <= j + (nKer-1)/2; n++) {
                    int c = n;
                    if (c < 0) c = 0;
                    if (c > nCol-1) c = nCol-1;
                    temp += (float)input[r][c]*kernel[kr][kc];
                    kc++;
                }
                kr++;
                kc = 0;
            }
            
            if (temp < 0) {
                temp = 0;}
            if (temp > 255) {
                temp = 255;}

            output[i][j] = (unsigned char)temp; 
        }
    }
}

void Image::BoxBlur(unsigned int kernel_size) {
    //create kernel matrix
    boost::multi_array<float, 2> kernel(boost::extents[kernel_size][kernel_size]);

    for (unsigned int i = 0; i < kernel_size; i++) {
        for (unsigned int j = 0; j < kernel_size; j++) {
            kernel[i][j] = (float)1.f/(kernel_size * kernel_size);
        }
    }
    Convolution(this->input, this->output, kernel);
}

unsigned int Image::Sharpness() {
    //create laplacian operator kernel
    boost::multi_array<float, 2> ker(boost::extents[3][3]);
    ker[0][0] = 0.f; ker[0][1] = 1.f; ker[0][2] = 0.f;
    ker[1][0] = 1.f; ker[1][1] = -4.f; ker[1][2] = 1.f;
    ker[2][0] = 0.f; ker[2][1] = 1.f; ker[2][2] = 0.f;

    boost::multi_array<unsigned char,2> result(boost::extents[this->output.shape()[0]][this->output.shape()[1]]);

    Convolution(this->output, result, ker);

    unsigned int sharpness = *std::max_element(result.data(), result.data() + result.num_elements());

    return sharpness;
}