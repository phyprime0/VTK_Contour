cmake_minimum_required(VERSION 3.14)

project(VTK_Contour LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew/")

# ---- Include VTK ----

find_package(VTK REQUIRED
    COMPONENTS
        CommonCore
        CommonDataModel
        FiltersExtraction
        FiltersGeometry
        FiltersSources
        FiltersAMR
        IOGeometry
        IOLegacy
        IOAMR
        RenderingCore
        TestingCore
        IOXML
)

if (NOT VTK_FOUND)
  message(FATAL_ERROR "VTK_Contour: Unable to find the VTK build folder.")
endif()



# ---- Declare library ----

add_library(
    VTK_Contour_lib
    source/vtk_contour.cpp
    include/vtk_contour.h
)
add_library(VTK_Contour::VTK_Contour ALIAS VTK_Contour_lib)

target_link_libraries(VTK_Contour_lib PRIVATE ${VTK_LIBRARIES})

vtk_module_autoinit(
    TARGETS VTK_Contour_lib
    MODULES ${VTK_LIBRARIES}
)

# include(GenerateExportHeader)
# generate_export_header(
#     build_test1_build_test1
#     BASE_NAME build_test1
#     EXPORT_FILE_NAME export/build_test1/build_test1_export.h
#     # CUSTOM_CONTENT_FROM_VARIABLE pragma_suppress_c4251
# )

if(NOT BUILD_SHARED_LIBS)
  target_compile_definitions(VTK_Contour_lib PUBLIC VTK_CONTOUR_STATIC_DEFINE)
endif()

set_target_properties(
    VTK_Contour_lib PROPERTIES
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN YES
    VERSION "${PROJECT_VERSION}"
    SOVERSION "${PROJECT_VERSION_MAJOR}"
    EXPORT_NAME VTK_Contour
    OUTPUT_NAME VTK_Contour
)

target_include_directories(
    VTK_Contour_lib ${warning_guard}
    PUBLIC
    "\$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
)

target_include_directories(
    VTK_Contour_lib SYSTEM
    PUBLIC
    "\$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/export>"
)

add_subdirectory(test)

target_compile_features(VTK_Contour_lib PUBLIC cxx_std_17)

# ---- Install rules ----

if(NOT CMAKE_SKIP_INSTALL_RULES)
  # include(cmake/install-rules.cmake)
endif()

# ---- Developer mode ----

if(NOT VTK_Contour_DEVELOPER_MODE)
  return()
elseif(NOT PROJECT_IS_TOP_LEVEL)
  message(
      AUTHOR_WARNING
      "Developer mode is intended for developers of build_test1"
  )
endif()

# include(cmake/dev-mode.cmake)


# add_executable(VTK_Contour main.cpp)

# include(GNUInstallDirs)
# install(TARGETS VTK_Contour
#     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
# )
