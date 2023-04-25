#ifndef ENGINE_MODEL
#define ENGINE_MODEL

#include "matrix.hpp"
#include "mesh.hpp"
#include "shader.hpp"

struct Model
{
    mat4 transform;

    int mesh_count;
    int material_count;

    Mesh *meshes;
    Material *materials;
    Material_Type *mesh_mateiral;
};

#endif