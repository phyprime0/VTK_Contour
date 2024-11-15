#include <string>
#include <iostream>
#include "vtk_contour.h"

auto main() -> int
{
  // auto const exported = exported_class {};

  // return std::string("build_test1") == exported.name() ? 0 : 1;


  auto *rtest = init_result();

  get_result(rtest,"test");

  size_t test_vertsize = get_verts_bytesize(rtest);
  size_t test_polysize = get_polys_bytesize(rtest);

  char* test_data = get_vertdata(rtest);

  float* f_data = reinterpret_cast<float*>(test_data);
  size_t num_v = test_vertsize / (sizeof(float) * 3);


  std::cout << "---------------------\n";
  std::cout << "Total verts: [" << num_v <<  "] \n Read (first 10) vert locations:\n";

  for (size_t i = 0; i < 10; i++){
      std::cout << "[" << f_data[i*3] << "," << f_data[1 + i*3] << "," << f_data[2 + i*3] << "]\n";
  }

  char* test_data2 = get_normdata(rtest);

  float* n_data = reinterpret_cast<float*>(test_data2);
  size_t num_n = test_vertsize / (sizeof(float) * 3);


  std::cout << "---------------------\n";
  std::cout << "Total norms: [" << num_n <<  "] \n Read (first 10) normal vectors:\n";

  for (size_t i = 0; i < 10; i++){
      std::cout << "[" << n_data[i*3] << "," << n_data[1 + i*3] << "," << n_data[2 + i*3] << "]\n";
  }

  char* test3_data = get_polydata(rtest);

  int* p_data = reinterpret_cast<int*>(test3_data);
  size_t num_p = test_polysize / (sizeof(int) * 3);

  std::cout << "---------------------\n";
  std::cout << "Total polys: [" << num_p <<  "] \n Read (first 10) poly indices:\n";

  for (size_t i = 0; i < 10; i++){
      std::cout << "[" << p_data[i*3] << "," << p_data[1 + i*3] << "," << p_data[2 + i*3] << "]\n";
  }


  deinit_result(rtest);

  return 0;

  // return 1;

}
