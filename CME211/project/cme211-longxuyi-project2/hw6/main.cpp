#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>
#include <jpeglib.h>
#include "hw6.hpp"
#include "image.hpp"

int main() {
  std::string filename = "stanford.jpeg";
  Image img(filename);

  std::cout << "Original Image: " << img.Sharpness() << std::endl;

  for (int i = 3; i <= 27; i= i+4) {
  Image img(filename);
  img.BoxBlur(i);
  std::cout << "Boxblur(" << i << "): "  << img.Sharpness() << std::endl;
  std::stringstream s;
  s << std::setw(2) << std::setfill('0') << i;
  std::string output_name = "BoxBlur" + s.str() + ".jpeg";
  img.Save(output_name);
  }
  
  return 0;
}