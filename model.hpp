#ifndef ENGINE_MODEL
#define ENGINE_MODEL

#include "mat.hpp"
#include "mesh.hpp"
#include "shader.hpp"

struct Model
{
    mat4 transform;

    int mesh_count;
    int material_count;

    Mesh *meshes;
    Shader *shaders;
    Material *materials;
    int *mesh_mateiral;
};

Model mesh_material_to_model(const Mesh &mesh, const Material &mateiral)
{
    Model model;

    model.transform = eye4();

    model.mesh_count = 1;
    model.meshes = (Mesh *)malloc(model.mesh_count * sizeof(Mesh));
    model.meshes[0] = mesh;

    model.shaders = (Shader *)malloc(model.mesh_count * sizeof(Shader));
    model.shaders[0].resize(model.meshes[0].triangle_count);

    model.material_count = 1;
    model.materials = (Material *)malloc(model.material_count * sizeof(Material));
    model.materials[0] = mateiral;

    model.mesh_mateiral = (int *)malloc(model.mesh_count * sizeof(int));
    model.mesh_mateiral[0] = 0;

    return model;
}

#endif