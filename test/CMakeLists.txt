cmake_minimum_required(VERSION 3.14)

project(build_test1Tests LANGUAGES CXX)

# include(../cmake/project-is-top-level.cmake)
# include(../cmake/folders.cmake)

# ---- Dependencies ----

if(PROJECT_IS_TOP_LEVEL)
  find_package(VTK_Contour REQUIRED)
  enable_testing()
endif()

# ---- Tests ----

add_executable(VTK_Contour_test source/vtk_contour_test.cpp)
target_link_libraries(VTK_Contour_test PRIVATE VTK_Contour::VTK_Contour)
target_compile_features(VTK_Contour_test PRIVATE cxx_std_17)

add_test(NAME build_test1_test COMMAND build_test1_test)

# ---- End-of-file commands ----

# add_folders(Test)
