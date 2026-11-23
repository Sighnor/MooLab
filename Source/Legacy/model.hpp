#ifndef MOOLAB_MODEL
#define MOOLAB_MODEL

#include <new>

#include "Core.h"
#include "mat.hpp"
#include "mesh.hpp"
#include "shader.hpp"

struct MooModel
{
    mat4 transform;

    int mesh_count;
    int material_count;

    MooMesh *meshes;
    MooShader *shaders;
    MooMaterial *materials;
    int *mesh_mateiral;
};

MooModel mesh_material_to_model(const MooMesh &mesh, const MooMaterial &mateiral)
{
    MooModel model;

    model.transform = eye4();

    model.mesh_count = 1;
    model.meshes = (MooMesh *)malloc(model.mesh_count * sizeof(MooMesh));
    model.meshes[0] = mesh;

    model.shaders = (MooShader *)malloc(model.mesh_count * sizeof(MooShader));
    MOO_CHECK_ALWAYS(model.shaders != NULL, "out of memory allocating shaders");
    //malloc只分配内存不构造对象，须先就地构造，否则resize读到的triangle_count/ver为未初始化值
    new (&model.shaders[0]) MooShader();
    model.shaders[0].resize(model.meshes[0].triangle_count);
    model.shaders[0].shader_type = mateiral.material_type;

    model.material_count = 1;
    model.materials = (MooMaterial *)malloc(model.material_count * sizeof(MooMaterial));
    model.materials[0] = mateiral;

    model.mesh_mateiral = (int *)malloc(model.mesh_count * sizeof(int));
    model.mesh_mateiral[0] = 0;

    return model;
}

#endif