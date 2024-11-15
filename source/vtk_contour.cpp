#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>



#include "vtk_contour.h"




#include <vtkAMReXGridReader.h>
#include <vtkOverlappingAMR.h>
#include <vtkUniformGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkCellDataToPointData.h>
#include <vtkDataSet.h>
#include <vtkAppendFilter.h>

#include <vtkDoubleArray.h>

#include <vtkContourFilter.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>

#include <vtkPolyDataWriter.h>


extern "C" {


struct Result{
  int num_verts;
  int num_polys;

  void* vert_data;
  void* poly_data;
  void* norm_data;

};

Result* init_result(){
  Result* result_data = new Result;

  result_data->num_polys = 0;
  result_data->num_verts = 0;

  return result_data;
};

void deinit_result(Result* res){
  delete res;
};


void get_result(Result* result_data, const char* filename) {
  // if (argc < 5) {
  //   return EXIT_FAILURE;
  // }

  // std::string input_fname = argv[1];
  // std::string var_name = argv[2];
  // double var_value = std::stod(argv[3]);
  // int fine_level = std::stoi(argv[4]);

  const std::string input_fname = "/Users/jwashin2/Documents/Projects/VTK/VTK-NOODLES/Data/plt00180";

  const std::string outputFileName = "/Users/jwashin2/Documents/Projects/VTK/VTK-NOODLES/Data/output.vtp";

  const std::string outputOBJName = "/Users/jwashin2/Documents/Projects/VTK/VTK-NOODLES/Data/output2.obj";

  const std::string var_name = "mag_vorticity";

  const int fine_level = 2;

  const double var_value = 0.2;

  auto reader = vtkSmartPointer<vtkAMReXGridReader>::New();

  reader->SetFileName(input_fname.c_str());
  reader->SetMaxLevel(fine_level);
  reader->Update(); // Without this, apparently it has no idea what arrays exist?? Does not match docs
  reader->SetCellArrayStatus(var_name.c_str(), 1);
  reader->UpdateInformation();
  reader->Update();

  auto amr_data = reader->GetOutput();

  std::cout << "AMR: points: " << amr_data->GetNumberOfPoints() << std::endl;
  std::cout << "AMR: cells: " << amr_data->GetNumberOfCells() << std::endl;
  std::cout << "AMR: levels: " << amr_data->GetNumberOfLevels() << std::endl;
  std::cout << "AMR: point arrays: " << reader->GetNumberOfPointArrays() << std::endl;
  std::cout << "AMR: cell arrays: " << reader->GetNumberOfCellArrays() << std::endl;

  std::cout << std::endl;
  auto amr_dataset = amr_data->GetDataSet(0,0);
  auto amr_cell_data = amr_dataset->GetCellData();

  assert(amr_cell_data->GetArray(var_name.c_str()));

  vtkSmartPointer<vtkAppendFilter> appendFilter = vtkSmartPointer<vtkAppendFilter>::New();

  vtkSmartPointer<vtkCompositeDataIterator> iter = amr_data->NewIterator();
  iter->InitTraversal();

  while (!iter->IsDoneWithTraversal()) {
    vtkDataSet* dataSet = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());
    if (dataSet) {
      appendFilter->AddInputData(dataSet);
    }
    iter->GoToNextItem();
  }

  appendFilter->Update();

  vtkSmartPointer<vtkUnstructuredGrid> ug_data = appendFilter->GetOutput();

  auto c2p = vtkSmartPointer<vtkCellDataToPointData>::New();
  c2p->SetInputData(ug_data);
  c2p->Update();
  auto c2p_result = c2p->GetOutput();

  assert(c2p_result);

  std::cout << "C2P: Points: " << c2p_result->GetNumberOfPoints() << std::endl;
  std::cout << "C2P: Cells: " << c2p_result->GetNumberOfCells() << std::endl;

  std::cout << std::endl;

  {
    auto pdata = c2p_result->GetPointData();
    assert(pdata->HasArray(var_name.c_str()));
    auto mvdata = pdata->GetArray(var_name.c_str());
    assert(mvdata);
    double range[2];
    mvdata->GetRange(range);
    std::cout << var_name << " range: " << range[0] << " - " << range[1] << std::endl;
  }

  auto contour = vtkSmartPointer<vtkContourFilter>::New();
  contour->SetInputData(c2p_result);
  contour->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, var_name.c_str());
  contour->SetValue(0, var_value);

  std::cout << "Contouring " << var_name << " at value: " << var_value << std::endl;

  std::cout << std::endl;

  //Contour filter will generate triangles and compute normals
  contour->SetGenerateTriangles(true);
  contour->SetComputeNormals(true);
  contour->Update();

  auto contour_result = contour->GetOutput();

  std::cout << "Contour: points: " << contour_result->GetNumberOfPoints() << std::endl;
  std::cout << "Contour: cells: " << contour_result->GetNumberOfCells() << std::endl;
  std::cout << "Contour: verts: " << contour_result->GetNumberOfVerts() << std::endl;
  std::cout << "Contour: polys: " << contour_result->GetNumberOfPolys() << std::endl;
  std::cout << "Contour: lines: " << contour_result->GetNumberOfLines() << std::endl;

  std::cout << std::endl;

  //Get all polygons as "cells" from the contour result
  vtkSmartPointer<vtkCellArray> polys = contour_result->GetPolys();

  auto cpoint_data = contour_result->GetPoints()->GetData();
  auto cvert_data = contour_result->GetVerts()->GetData();


  auto cnorms = contour_result->GetPointData()->GetNormals();

  //Get all polygon points from contour result

  auto num_vals = cpoint_data->GetNumberOfValues();
  auto num_tups = cpoint_data->GetNumberOfTuples();
  auto num_comps = cpoint_data->GetNumberOfComponents();

  auto num_ntups = cnorms->GetNumberOfTuples();


  //Set up 2D vector for point coordinates
  std::vector<std::array<float,3>> poly_coords(num_tups);

  std::vector<std::array<float,3>> poly_norms(num_ntups);

  //Iterate through points and store coordinates/normals in vector
  for (auto tuple_id = 0; tuple_id < num_tups; tuple_id++){
    for (auto c = 0; c < num_comps; c++){

        poly_coords[tuple_id][c] = cpoint_data->GetComponent(tuple_id,c);
        poly_norms[tuple_id][c] = cnorms->GetComponent(tuple_id,c);
    }
  }

  //Set up iteration for polygons (as "cells")
  vtkIdType cellID = 1;
  //IdList serves to store reference ID for points of polygon
  auto idList = vtkSmartPointer<vtkIdList>::New();
  auto count = 0;
  int total_cells = polys->GetNumberOfCells();

  //Set up 2D vector for polygons (vertex indices), made variable for polygons > tris, but currently coded for tris only.
  std::vector<std::array<int,3>> poly_indices(total_cells);

  //Iterate through polygon "cells" to extract vertex indices into vectors
  for (vtkIdType i = 0; i < polys->GetNumberOfCells(); i++) {
    polys->GetCellAtId(i,idList);
    for (auto k = 0; k < idList->GetNumberOfIds(); k++) {
      poly_indices[i][k] = idList->GetId(k);
    }
    if (idList->GetNumberOfIds()>3) {
      std::cout << "non-tri: " << idList->GetNumberOfIds() << std::endl;
    }

  }


  if (true){
    //Write to .OBJ file for debugging/testomg
    std::ofstream ofile(outputOBJName);
    ofile << "#OBJ File" << "\n" << "o VTK_Contour" << "\n";
    //Write all vertices
    for (auto l = 0; l < poly_coords.size(); l++) {
      ofile << "v " << poly_coords[l][0] << " " << poly_coords[l][1] << " " << poly_coords[l][2] << "\n";
    }
    //Write all faces
    ofile << "\n";
    for (auto l = 0; l < poly_indices.size(); l++) {
      ofile << "f " << poly_indices[l][0]+1 << " " << poly_indices[l][1]+1 << " " << poly_indices[l][2]+1 << "\n";
    }
    ofile.close();
    std::cout << "File writing complete" << std::endl;
  }

  result_data->num_polys = total_cells;
  result_data->num_verts = num_tups;


  //VERT DATA

  size_t byte_count1 = num_tups * 3 * sizeof(float);

  char* ptr1 = new char[byte_count1];

  size_t index1 = 0;

  for (const auto& row : poly_coords) {
      for (float value : row) {
          std::memcpy(ptr1 + index1, &value, sizeof(float));
          index1 += sizeof(float);  // Move the index forward by the size of one integer
      }
  }


  //POLY DATA

  size_t byte_count2 = total_cells * 3 * sizeof(int);

  char* ptr2 = new char[byte_count2];

  size_t index2 = 0;

  for (const auto& row : poly_indices) {
      for (int value : row) {
          std::memcpy(ptr2 + index2, &value, sizeof(int));
          index2 += sizeof(int);  // Move the index forward by the size of one integer
      }
  }


  //NORMALS DATA

  size_t byte_count3 = num_ntups * 3 * sizeof(float);

  char* ptr3 = new char[byte_count3];

  size_t index3 = 0;

  for (const auto& row : poly_norms) {
      for (float value : row) {
          std::memcpy(ptr3 + index3, &value, sizeof(float));
          index3 += sizeof(float);  // Move the index forward by the size of one integer
      }
  }


  result_data->vert_data = ptr1;

  result_data->poly_data = ptr2;

  result_data->norm_data = ptr3;

  // result_data->poly_data = poly_indices.data();
  // result_data->vert_data = poly_coords.data();



  // result_data->poly_data = poly_indices;
  // result_data->vert_data = poly_coords;


}




// int get_num_verts(Result* res){
//   return res->num_verts;
// }

// int get_num_polys(Result* res){
//   return res->num_polys;
// }

size_t get_polys_bytesize(Result* res) {
  size_t s_result = res->num_polys;
  size_t byte_count = s_result * 3 * sizeof(int);

  std::cout << "Polys: " << s_result << " [Bytes]: " << byte_count << std::endl;

  return byte_count;
}

size_t get_verts_bytesize(Result* res) {
  size_t s_result = res->num_verts;
  size_t byte_count = s_result * 3 * sizeof(float);

  std::cout << "Verts: " << s_result << " [Bytes]: " << byte_count << std::endl;

  return byte_count;
}

char* get_vertdata(Result* res) {
    // size_t v_size = get_verts_bytesize(res);

    char* ptr = static_cast<char*>(res->vert_data);
    // std::memcpy(ptr1, res->vert_data, v_size);

    return ptr;
}


char* get_normdata(Result* res) {
    // size_t v_size = get_verts_bytesize(res);

    char* ptr = static_cast<char*>(res->norm_data);
    // std::memcpy(ptr1, res->vert_data, v_size);
    return ptr;
}


char* get_polydata(Result* res) {
    // size_t p_size = get_polys_bytesize(res);

    // char* ptr2 =new char[p_size];
    // std::memcpy(ptr2, res->poly_data, p_size);

    // return ptr2;

    char* ptr = static_cast<char*>(res->poly_data);
    return ptr;
}


}

