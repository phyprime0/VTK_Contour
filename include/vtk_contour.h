#pragma once

#include <cstdlib>

// extern "C" {

struct Configuration;

Configuration* new_configuration();
void           destroy_configuration(Configuration*);

void set_input_filename(Configuration*, const char* file);
void set_variable_name(Configuration*, const char* name);
void set_variable_value(Configuration*, double value);
void set_max_refinement_level(Configuration*, unsigned level);


struct Result;

Result* init_result(Configuration*);
void    deinit_result(Result*);

char* get_vertdata(Result*);
char* get_polydata(Result*);
char* get_normdata(Result*);

size_t get_polys_bytesize(Result*);
size_t get_verts_bytesize(Result*);

void dump_to_obj(Result*, const char* file);
// }
