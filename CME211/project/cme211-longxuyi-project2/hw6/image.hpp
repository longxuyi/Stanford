#ifndef IMAGE_HPP
#define IMAGE_HPP

#include <string>
#include <boost/multi_array.hpp>

class Image {

    std::string input_file;
    boost::multi_array<unsigned char, 2> input;
    boost::multi_array<unsigned char, 2> output;
    boost::multi_array<unsigned char, 2> img;


 public:

    //Image() constructor.
    Image(std::string filename);

    // Save() writes current version of the image to a jpeg file.
    void Save(std::string filename);

    //the input and output should be of the same size and only support odd size
    void Convolution(boost::multi_array<unsigned char,2>& input, 
                 boost::multi_array<unsigned char,2>& output,
                 boost::multi_array<float,2>& kernel);

    /* Boxblur() 
    input arg: the kernel size, output: uses Convolution to smooth the image. The values in the kernel 
    should always be scaled so that they sum to 1 */
    void BoxBlur(unsigned int kernel_size);


    /* Sharpness()
    returns unsigned int for the sharpness of the image. Use equation 2 and Convolution()
    function to return the max of convolution output. */
    unsigned int Sharpness();

};

#endif /* IMAGE_HPP */