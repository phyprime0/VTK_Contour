// #include <string>
#include <iostream>
#include "vtk_contour.h"

auto main() -> int {

    auto* config = new_configuration();

    set_input_filename(
        config,
        "/Users/jwashin2/Documents/Projects/VTK/VTK-NOODLES/Data/plt00180");
    set_max_refinement_level(config, 2);
    set_variable_name(config, "mag_vorticity");
    set_variable_value(config, 0.2);
    set_texture_variable(config, "p");

    auto* rtest = init_result(config);

    destroy_configuration(config);

    size_t test_vertsize = get_verts_bytesize(rtest);
    size_t test_polysize = get_polys_bytesize(rtest);
    size_t test_texturesize = get_texture_bytesize(rtest);

    char* test_data = get_vertdata(rtest);

    float* f_data = reinterpret_cast<float*>(test_data);
    size_t num_v  = test_vertsize / (sizeof(float) * 3);


    std::cout << "---------------------\n";
    std::cout << "Total verts: [" << num_v
              << "] \n Read (first 10) vert locations:\n";

    for (size_t i = 0; i < 10; i++) {
        std::cout << "[" << f_data[i * 3] << "," << f_data[1 + i * 3] << ","
                  << f_data[2 + i * 3] << "]\n";
    }

    char* test_data2 = get_normdata(rtest);

    float* n_data = reinterpret_cast<float*>(test_data2);
    size_t num_n  = test_vertsize / (sizeof(float) * 3);


    std::cout << "---------------------\n";
    std::cout << "Total norms: [" << num_n
              << "] \n Read (first 10) normal vectors:\n";

    for (size_t i = 0; i < 10; i++) {
        std::cout << "[" << n_data[i * 3] << "," << n_data[1 + i * 3] << ","
                  << n_data[2 + i * 3] << "]\n";
    }

    char* test3_data = get_polydata(rtest);

    int*   p_data = reinterpret_cast<int*>(test3_data);
    size_t num_p  = test_polysize / (sizeof(int) * 3);

    std::cout << "---------------------\n";
    std::cout << "Total polys: [" << num_p
              << "] \n Read (first 10) poly indices:\n";

    for (size_t i = 0; i < 10; i++) {
        std::cout << "[" << p_data[i * 3] << "," << p_data[1 + i * 3] << ","
                  << p_data[2 + i * 3] << "]\n";
    }


    char* test4_data = get_texcoorddata(rtest);

    float*   t_data = reinterpret_cast<float*>(test4_data);
    size_t num_t  = test_polysize / (sizeof(float) * 1);

    std::cout << "---------------------\n";
    std::cout << "Total Texture_Coordinates: [" << num_t
              << "] \n Read (first 10) texture coordinates:\n";

    for (size_t i = 0; i < 10; i++) {
        std::cout << "[" << t_data[i] << "]\n";
    }









    deinit_result(rtest);

    return EXIT_SUCCESS;
}
