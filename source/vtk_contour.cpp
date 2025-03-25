#include <cassert>
#include <cstdlib>
#include <cstring> // for memcpy
#include <fstream>
#include <iostream>
#include <string>

#include "vtk_contour.h"

#include <vtkAMReXGridReader.h>
#include <vtkAppendFilter.h>
#include <vtkCellDataToPointData.h>
#include <vtkDataSet.h>
#include <vtkOverlappingAMR.h>
#include <vtkSmartPointer.h>
#include <vtkUniformGrid.h>
#include <vtkUnstructuredGrid.h>

#include <vtkDoubleArray.h>
#include <vtkColorSeries.h>
#include <vtkImageData.h>
#include <vtkTexture.h>
#include <vtkDataArray.h>
#include <vtkPolyDataMapper.h>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkContourFilter.h>
#include <vtkPointData.h>

#include <vtkProbeFilter.h>

#include <vtkPolyDataWriter.h>


// extern "C" {

struct Configuration {
    std::string input_filename;
    std::string variable_name;
    std::string texture_variable;
    double      variable_value;
    unsigned    max_refinements;
};

Configuration* new_configuration() {
    return new Configuration();
}

void destroy_configuration(Configuration* c) {
    delete c;
}

void set_input_filename(Configuration* c, const char* file) {
    if (!c) return;
    c->input_filename = file;
}

void set_variable_name(Configuration* c, const char* name) {
    if (!c) return;
    c->variable_name = name;
}

void set_variable_value(Configuration* c, double value) {
    if (!c) return;
    c->variable_value = value;
}

void set_texture_variable(Configuration* c, const char* tname) {
    if (!c) return;
    c->texture_variable = tname;
}

void set_max_refinement_level(Configuration* c, unsigned level) {
    if (!c) return;
    c->max_refinements = level;
}


// =============================================================================

struct Result {
    int num_verts = 0;
    int num_polys = 0;
    //Basic colormap size, can be adjusted if needed.
    int texture_size = 256;


    std::vector<std::array<float, 3>>    vert_data;
    std::vector<std::array<float, 3>>    norm_data;
    std::vector<std::array<uint32_t, 3>> poly_data;
    std::vector<std::array<uint8_t, 3>>  texture_data;
    std::vector<float>                   texcoord_data;
};

Result* init_result(Configuration* configuration) {
    Result* result_data = new Result;

    const std::string input_fname = configuration->input_filename;
    const std::string var_name    = configuration->variable_name;
    const int         fine_level  = configuration->max_refinements;
    const double      var_value   = configuration->variable_value;
    const std::string tex_var     = configuration->texture_variable;


    auto reader = vtkSmartPointer<vtkAMReXGridReader>::New();

    reader->SetFileName(input_fname.c_str());
    reader->SetMaxLevel(fine_level);
    reader->Update(); // Without this, apparently it has no idea what arrays
                      // exist?? Does not match docs
    reader->SetCellArrayStatus(var_name.c_str(), 1);
    // Retrieve the texture variable by name and set ArrayStatus to True
    reader->SetCellArrayStatus(tex_var.c_str(), 1);
    reader->UpdateInformation();
    reader->Update();

    auto amr_data = reader->GetOutput();

    std::cout << "AMR: points: " << amr_data->GetNumberOfPoints() << std::endl;
    std::cout << "AMR: cells: " << amr_data->GetNumberOfCells() << std::endl;
    std::cout << "AMR: levels: " << amr_data->GetNumberOfLevels() << std::endl;
    std::cout << "AMR: point arrays: " << reader->GetNumberOfPointArrays()
              << std::endl;
    std::cout << "AMR: cell arrays: " << reader->GetNumberOfCellArrays()
              << std::endl;

    std::cout << std::endl;
    auto amr_dataset   = amr_data->GetDataSet(0, 0);
    auto amr_cell_data = amr_dataset->GetCellData();

    assert(amr_cell_data->GetArray(var_name.c_str()));
    assert(amr_cell_data->GetArray(tex_var.c_str()));

    auto appendFilter = vtkSmartPointer<vtkAppendFilter>::New();

    vtkSmartPointer<vtkCompositeDataIterator> iter = amr_data->NewIterator();
    iter->InitTraversal();

    while (!iter->IsDoneWithTraversal()) {
        vtkDataSet* dataSet =
            vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());
        if (dataSet) { appendFilter->AddInputData(dataSet); }
        iter->GoToNextItem();
    }

    appendFilter->Update();

    vtkSmartPointer<vtkUnstructuredGrid> ug_data = appendFilter->GetOutput();

    auto c2p = vtkSmartPointer<vtkCellDataToPointData>::New();
    c2p->SetInputData(ug_data);
    c2p->Update();
    auto c2p_result = c2p->GetOutput();

    assert(c2p_result);

    std::cout << "C2P: Points: " << c2p_result->GetNumberOfPoints()
              << std::endl;
    std::cout << "C2P: Cells: " << c2p_result->GetNumberOfCells() << std::endl;

    std::cout << std::endl;

    double txrange[2];
    // auto txdata = vtkSmartPointer<vtkDataArray>::New();

    vtkSmartPointer<vtkDataArray> txdata;

    {
        auto pdata = c2p_result->GetPointData();
        assert(pdata->HasArray(var_name.c_str()));
        auto mvdata = pdata->GetArray(var_name.c_str());
        assert(mvdata);
        double range[2];
        mvdata->GetRange(range);
        std::cout << var_name << " range: " << range[0] << " - " << range[1]
                  << std::endl;

        assert(pdata->HasArray(tex_var.c_str()));
        // auto txdata = pdata->GetArray(tex_var.c_str());
        txdata = pdata->GetArray(tex_var.c_str());
        assert(txdata);

        txdata->GetRange(txrange);
        std::cout << tex_var << " range: " << txrange[0] << " - " << txrange[1]
                  << std::endl;
    }


    // Generate Colormap (Might be smarter to restructure this... but for now...

    vtkSmartPointer<vtkColorSeries> colorSeries = vtkSmartPointer<vtkColorSeries>::New();
    //Can switch color schemes
    colorSeries->SetColorScheme(vtkColorSeries::BREWER_DIVERGING_PURPLE_ORANGE_11);

    int num_colors = colorSeries->GetNumberOfColors();

    int tx_size = result_data->texture_size;

    std::vector<std::array<uint8_t, 3>> c_texture(tx_size);


    // Interpolate colors manually
    for (int i = 0; i < tx_size; i++) {

        uint8_t r;
        uint8_t g;
        uint8_t b;


        double c = i / double(tx_size - 1);

        double c_index = double(c * (num_colors - 1));

        int intColor1 = std::floor(c_index);
        int intColor2 = intColor1+1;

        if (i == 0) {
            vtkColor3ub color = colorSeries->GetColor(0);
            r = color.GetRed();
            g = color.GetGreen();
            b = color.GetBlue();


        } else if (i == tx_size-1) {
            vtkColor3ub color = colorSeries->GetColor(num_colors-1);

            r = color.GetRed();
            g = color.GetGreen();
            b = color.GetBlue();

        } else {

            double c_frac = c_index - intColor1;
            double c_comp = 1 - c_frac;

            vtkColor3ub color1 = colorSeries->GetColor(intColor1);
            vtkColor3ub color2 = colorSeries->GetColor(intColor2);

            r = (color1.GetRed() * c_comp + color2.GetRed() * c_frac);
            g = (color1.GetGreen() * c_comp + color2.GetGreen() * c_frac);
            b = (color1.GetBlue() * c_comp + color2.GetBlue() * c_frac);

        }

        c_texture[i][0] = r;
        c_texture[i][1] = g;
        c_texture[i][2] = b;

    }

    result_data->texture_data = std::move(c_texture);





    auto contour = vtkSmartPointer<vtkContourFilter>::New();
    contour->SetInputData(c2p_result);
    contour->SetInputArrayToProcess(
        0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, var_name.c_str());
    contour->SetValue(0, var_value);

    std::cout << "Contouring " << var_name << " at value: " << var_value
              << std::endl;

    std::cout << std::endl;

    // Contour filter will generate triangles and compute normals
    contour->SetGenerateTriangles(true);
    contour->SetComputeNormals(true);
    contour->Update();

    auto contour_result = contour->GetOutput();

    std::cout << "Contour: points: " << contour_result->GetNumberOfPoints()
              << std::endl;
    std::cout << "Contour: cells: " << contour_result->GetNumberOfCells()
              << std::endl;
    std::cout << "Contour: verts: " << contour_result->GetNumberOfVerts()
              << std::endl;
    std::cout << "Contour: polys: " << contour_result->GetNumberOfPolys()
              << std::endl;
    std::cout << "Contour: lines: " << contour_result->GetNumberOfLines()
              << std::endl;

    std::cout << std::endl;


    // Probe filter will map the second variable onto the vertices of the countour

    vtkSmartPointer<vtkProbeFilter> probe = vtkSmartPointer<vtkProbeFilter>::New();
    probe->SetSourceData(c2p_result);  // The original dataset with the second cell array
    probe->SetInputData(contour_result);
    probe->Update();

    auto probe_result = probe->GetOutput();

    // vtkSmartPointer<vtkPolyDataMapper> probe = vtkSmartPointer<vtkPolyDataMapper>::New();
    // probe->  // The original dataset with the second cell array
    // probe->SetInputData(contour_result);
    // probe->Update();





    // Get all texture-variable values
    auto tpoint_data = probe_result->GetPointData();


    std::vector<float> poly_tdata(probe_result->GetNumberOfPoints());

    double abs_trange = abs(txrange[1]) + abs(txrange[0]);


    for (vtkIdType tid = 0; tid < probe_result->GetNumberOfPoints(); tid++) {
        auto tx_array = tpoint_data->GetArray(tex_var.c_str());
        auto tx_val = tx_array->GetComponent(tid, 0);

        auto tx_scaled = (tx_val- txrange[0])/(abs_trange);

        poly_tdata[tid] = tx_scaled;

    }





    // Get all polygons as "cells" from the contour result
    vtkSmartPointer<vtkCellArray> polys = contour_result->GetPolys();

    auto cpoint_data = contour_result->GetPoints()->GetData();
    auto cvert_data  = contour_result->GetVerts()->GetData();


    auto cnorms = contour_result->GetPointData()->GetNormals();


    // Get all polygon points from contour result

    auto num_vals  = cpoint_data->GetNumberOfValues();
    auto num_tups  = cpoint_data->GetNumberOfTuples();
    auto num_comps = cpoint_data->GetNumberOfComponents();

    auto num_ntups = cnorms->GetNumberOfTuples();


    // Set up 2D vector for point coordinates
    std::vector<std::array<float, 3>> poly_coords(num_tups);

    std::vector<std::array<float, 3>> poly_norms(num_ntups);

    // Iterate through points and store coordinates/normals in vector
    for (auto tuple_id = 0; tuple_id < num_tups; tuple_id++) {
        for (auto c = 0; c < num_comps; c++) {

            poly_coords[tuple_id][c] = cpoint_data->GetComponent(tuple_id, c);
            poly_norms[tuple_id][c]  = cnorms->GetComponent(tuple_id, c);
        }
    }

    // Set up iteration for polygons (as "cells")
    vtkIdType cellID = 1;
    // IdList serves to store reference ID for points of polygon
    auto idList      = vtkSmartPointer<vtkIdList>::New();
    auto count       = 0;
    int  total_cells = polys->GetNumberOfCells();

    // Set up 2D vector for polygons (vertex indices), made variable for
    // polygons > tris, but currently coded for tris only.
    std::vector<std::array<uint32_t, 3>> poly_indices(total_cells);

    // Iterate through polygon "cells" to extract vertex indices into vectors
    for (vtkIdType i = 0; i < polys->GetNumberOfCells(); i++) {
        polys->GetCellAtId(i, idList);
        for (auto k = 0; k < idList->GetNumberOfIds(); k++) {
            poly_indices[i][k] = idList->GetId(k);
        }
        if (idList->GetNumberOfIds() > 3) {
            std::cout << "non-tri: " << idList->GetNumberOfIds() << std::endl;
        }
    }

    result_data->num_polys = total_cells;
    result_data->num_verts = num_tups;

    result_data->vert_data = std::move(poly_coords);

    result_data->norm_data = std::move(poly_norms);

    result_data->poly_data = std::move(poly_indices);

    result_data->texcoord_data = std::move(poly_tdata);


    // result_data->poly_data = poly_indices.data();
    // result_data->vert_data = poly_coords.data();


    // result_data->poly_data = poly_indices;
    // result_data->vert_data = poly_coords;

    return result_data;
};

void deinit_result(Result* res) {
    delete res;
};

// int get_num_verts(Result* res){
//   return res->num_verts;
// }

// int get_num_polys(Result* res){
//   return res->num_polys;
// }

size_t get_polys_bytesize(Result* res) {
    size_t s_result   = res->num_polys;
    size_t byte_count = s_result * 3 * sizeof(int);

    std::cout << "Polys: " << s_result << " [Bytes]: " << byte_count
              << std::endl;

    return byte_count;
}

size_t get_verts_bytesize(Result* res) {
    size_t s_result   = res->num_verts;
    size_t byte_count = s_result * 3 * sizeof(float);

    std::cout << "Verts: " << s_result << " [Bytes]: " << byte_count
              << std::endl;

    return byte_count;
}

size_t get_texture_bytesize(Result* res) {
    size_t s_result   = res->texture_size;
    size_t byte_count = s_result * 3 * sizeof(uint8_t);

    std::cout << "Texture Data Size: " << s_result << " [Bytes]: " << byte_count
              << std::endl;

    return byte_count;
}

size_t get_texcoord_bytesize(Result* res) {
    size_t s_result   = res->num_verts;
    size_t byte_count = s_result * sizeof(float);

    std::cout << "Texture Coordinates: " << s_result << " [Bytes]: " << byte_count
              << std::endl;

    return byte_count;
}


char* get_vertdata(Result* res) {
    // size_t v_size = get_verts_bytesize(res);

    char* ptr = reinterpret_cast<char*>(res->vert_data.data());
    // std::memcpy(ptr1, res->vert_data, v_size);

    return ptr;
}


char* get_normdata(Result* res) {
    // size_t v_size = get_verts_bytesize(res);

    char* ptr = reinterpret_cast<char*>(res->norm_data.data());
    // std::memcpy(ptr1, res->vert_data, v_size);
    return ptr;
}


char* get_polydata(Result* res) {
    // size_t p_size = get_polys_bytesize(res);

    // char* ptr2 =new char[p_size];
    // std::memcpy(ptr2, res->poly_data, p_size);

    // return ptr2;

    char* ptr = reinterpret_cast<char*>(res->poly_data.data());
    return ptr;
}

char* get_texturedata(Result* res) {
    char* ptr = reinterpret_cast<char*>(res->texture_data.data());
    return ptr;
}

char* get_texcoorddata(Result* res) {
    char* ptr = reinterpret_cast<char*>(res->texcoord_data.data());
    return ptr;
}



void dump_to_obj(Result* res, const char* file) {
    // Write to .OBJ file for debugging/testomg
    std::ofstream ofile(file);
    ofile << "#OBJ File" << "\n" << "o VTK_Contour" << "\n";

    auto const& coords  = res->vert_data;
    auto const& normals = res->norm_data;
    auto const& index   = res->poly_data;

    // Write all vertices
    for (auto const& coord : coords) {
        ofile << "v " << coord[0] << " " << coord[1] << " " << coord[2] << "\n";
    }

    ofile << "\n";

    // Write all normals
    for (auto const& norm : normals) {
        ofile << "v " << norm[0] << " " << norm[1] << " " << norm[2] << "\n";
    }

    ofile << "\n";

    // Write all faces
    for (auto const& tri : index) {
        ofile << "f " << tri[0] + 1 << " " << tri[1] + 1 << " " << tri[2] + 1
              << "\n";
    }
    ofile.close();
}
// }
