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

    float vertices[9] = {
            -1.f, 0.f, 0.f,
            1.f, 0.f, 0.f,
            0.f, 1.f, 0.f};

    float normals[9] = {
            -1.f, 0.f, 0.f,
            1.f, 0.f, 0.f,
            0.f, 1.f, 0.f};

    float texcoords[6] = {
            0.f, 0.f,
            0.f, 0.f,
            0.f, 0.f};

    unsigned short indices[3] = {
            0, 1, 2};

    mesh.vertex_count = 3;
    mesh.triangle_count = 1;
    mesh.vertices = (float *)malloc(3 * 3 * sizeof(float));
    mesh.normals = (float *)malloc(3 * 3 * sizeof(float));
    mesh.texcoords = (float *)malloc(3 * 2 * sizeof(float));
    mesh.indices = (unsigned short *)malloc(1 * 3 * sizeof(unsigned short));

    memcpy(mesh.vertices, vertices, 3 * 3 * sizeof(float));
    memcpy(mesh.normals, normals, 3 * 3 * sizeof(float));
    memcpy(mesh.texcoords, texcoords, 3 * 2 * sizeof(float));
    memcpy(mesh.indices, indices, 1 * 3 * sizeof(unsigned short));

    return mesh;
}

Mesh light_cube(mat3 orientation, vec3 position)
{
    Mesh mesh;

    float vertices[72] = {
            // Front face
            -1.f, -1.f, 1.f,
            1.f, -1.f, 1.f,
            1.f, 1.f, 1.f,
            -1.f, 1.f, 1.f,

            // Back face
            -1.f, -1.f, -1.f,
            -1.f, 1.f, -1.f,
            1.f, 1.f, -1.f,
            1.f, -1.f, -1.f,

            // Top face
            -1.f, 1.f, -1.f,
            -1.f, 1.f, 1.f,
            1.f, 1.f, 1.f,
            1.f, 1.f, -1.f,

            // Bottom face
            -1.f, -1.f, -1.f,
            1.f, -1.f, -1.f,
            1.f, -1.f, 1.f,
            -1.f, -1.f, 1.f,

            // Right face
            1.f, -1.f, -1.f,
            1.f, 1.f, -1.f,
            1.f, 1.f, 1.f,
            1.f, -1.f, 1.f,

            // Left face
            -1.f, -1.f, -1.f,
            -1.f, -1.f, 1.f,
            -1.f, 1.f, 1.f,
            -1.f, 1.f, -1.f};

    float normals[72] = {
            // Front face
            0.f, 0.f, 1.f,
            0.f, 0.f, 1.f,
            0.f, 0.f, 1.f,
            0.f, 0.f, 1.f,

            // Back face
            0.f, 0.f, -1.f,
            0.f, 0.f, -1.f,
            0.f, 0.f, -1.f,
            0.f, 0.f, -1.f,

            // Top face
            0.f, 1.f, 0.f,
            0.f, 1.f, 0.f,
            0.f, 1.f, 0.f,
            0.f, 1.f, 0.f,

            // Bottom face
            0.f, -1.f, 0.f,
            0.f, -1.f, 0.f,
            0.f, -1.f, 0.f,
            0.f, -1.f, 0.f,

            // Right face
            1.f, 0.f, 0.f,
            1.f, 0.f, 0.f,
            1.f, 0.f, 0.f,
            1.f, 0.f, 0.f,

            // Left face
            -1.f, 0.f, 0.f,
            -1.f, 0.f, 0.f,
            -1.f, 0.f, 0.f,
            -1.f, 0.f, 0.f};

    float texcoords[48] = {
            // Front face
            0.f, 0.f,
            1.f, 0.f,
            1.f, 1.f,
            0.f, 1.f,

            // Back face
            0.f, 0.f,
            1.f, 0.f,
            1.f, 1.f,
            0.f, 1.f,

            // Top face
            0.f, 0.f,
            1.f, 0.f,
            1.f, 1.f,
            0.f, 1.f,

            // Bottom face
            0.f, 0.f,
            1.f, 0.f,
            1.f, 1.f,
            0.f, 1.f,

            // Right face
            0.f, 0.f,
            1.f, 0.f,
            1.f, 1.f,
            0.f, 1.f,

            // Left face
            0.f, 0.f,
            1.f, 0.f,
            1.f, 1.f,
            0.f, 1.f};

    unsigned short indices[36] = {
            0, 1, 2, 0, 2, 3,
            4, 5, 6, 4, 6, 7,
            8, 9, 10, 8, 10, 11,
            12, 13, 14, 12, 14, 15,
            16, 17, 18, 16, 18, 19,
            20, 21, 22, 20, 22, 23};

    mesh.vertex_count = 24;
    mesh.triangle_count = 12;
    mesh.vertices = (float *)malloc(24 * 3 * sizeof(float));
    mesh.normals = (float *)malloc(24 * 3 * sizeof(float));
    mesh.texcoords = (float *)malloc(24 * 2 * sizeof(float));
    mesh.indices = (unsigned short *)malloc(12 * 3 * sizeof(unsigned short));

    memcpy(mesh.vertices, vertices, 24 * 3 * sizeof(float));
    memcpy(mesh.normals, normals, 24 * 3 * sizeof(float));
    memcpy(mesh.texcoords, texcoords, 24 * 2 * sizeof(float));
    memcpy(mesh.indices, indices, 12 * 3 * sizeof(unsigned short));

    return mesh;
}

#endif