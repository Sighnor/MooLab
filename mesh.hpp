#ifndef ENGINE_MESH
#define ENGINE_MESH

#include "mat.hpp"

struct Mesh
{
    int vertex_count;
    int triangle_count;
    float* vertices;
    float* normals;
    float* texcoords;
    unsigned short* indices;
};

Mesh triangle(mat3 orientation, vec3 position)
{
    Mesh mesh;

    float vertices[9] = {-1.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 1.f, 0.f};
    float normals[9] = {-1.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 1.f, 0.f};
    float texcoords[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
    unsigned short indices[9] = {0, 1, 2};

    mesh.vertex_count = 3;
    mesh.triangle_count = 1;
    mesh.vertices = (float *)malloc(3 * 3 * sizeof(float));
    mesh.normals = (float *)malloc(3 * 3 * sizeof(float));
    mesh.texcoords = (float *)malloc(3 * 2 * sizeof(float));
    mesh.indices = (unsigned short *)malloc(3 * sizeof(unsigned short));

    memcpy(mesh.vertices, vertices, 3 * 3 * sizeof(float));
    memcpy(mesh.normals, normals, 3 * 3 * sizeof(float));
    memcpy(mesh.texcoords, texcoords, 3 * 2 * sizeof(float));
    memcpy(mesh.indices, indices, 3 * sizeof(unsigned short));

    return mesh;
}

Mesh light_cube(mat3 orientation, vec3 position)
{
    Mesh mesh;

    mesh.vertex_count = 8;
    mesh.triangle_count = 12;
    mesh.vertices = (float *)malloc(8 * 3 * sizeof(float));
    mesh.normals = (float *)malloc(12 * 3 * sizeof(float));
    mesh.texcoords = (float *)malloc(12 * 2 * sizeof(float));
    mesh.indices = (unsigned short *)malloc(12 * sizeof(unsigned short));

    return mesh;
}

#endif