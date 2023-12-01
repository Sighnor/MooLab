#ifndef MOOLAB_MESH
#define MOOLAB_MESH

#include "mat.hpp"

struct MooMesh
{
    int vertex_count;
    int triangle_count;
    float* vertices;
    float* normals;
    float* texcoords;
    unsigned short* indices;
};

MooMesh triangle(mat3 orientation, vec3 position)
{
    MooMesh mesh;

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

MooMesh cube()
{
    MooMesh mesh;

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

MooMesh capsule(float r, float h, int r_size, int h_size, float s = 0)
{
        MooMesh mesh;

        int triangles_size = 2 * r_size + 2 * h_size * r_size;
        int vertices_size = 2 + (h_size + 1) * r_size;

        float vertices[3 * vertices_size];

        float normals[3 * vertices_size];

        float texcoords[2 * vertices_size];

        unsigned short indices[3 * triangles_size];

        vertices[0] = 0.f;
        vertices[1] = -h / 2.f;
        vertices[2] = 0.f;
        normals[0] = 0.f;
        normals[1] = -1.f;
        normals[2] = 0.f;
        texcoords[0] = 0.f;
        texcoords[1] = 0.f;

        vertices[3 * (vertices_size - 1)] = 0.f;
        vertices[3 * (vertices_size - 1) + 1] = h / 2.f;
        vertices[3 * (vertices_size - 1) + 2] = 0.f;
        normals[3 * (vertices_size - 1)] = 0.f;
        normals[3 * (vertices_size - 1) + 1] = 1.f;
        normals[3 * (vertices_size - 1) + 2] = 0.f;
        texcoords[2 * (vertices_size - 1)] = 0.f;
        texcoords[2 * (vertices_size - 1) + 1] = 0.f;

        for(int i = 0; i < (h_size + 1); i++)
        {
                for(int j = 0; j < r_size; j++)
                {
                        int id = 1 + i * r_size + j;
                        float phi = j * 2.f * PI / r_size;
                        float k = float(i) / h_size;
                        vertices[3 * id] = r * sin(phi);
                        vertices[3 * id + 1] = h * k - h / 2.f;
                        vertices[3 * id + 2] = r * cos(phi);
                        if(abs(vertices[3 * id + 1]) > h / 2.f - s)
                        {
                                float scale = sin(acos((abs(vertices[3 * id + 1]) - h / 2.f + s) / std::max(s, 0.001f)));
                                vertices[3 * id] = scale * vertices[3 * id];
                                vertices[3 * id + 2] = scale * vertices[3 * id + 2];
                        }
                        normals[3 * id] = sin(phi);
                        normals[3 * id + 1] = 0.f;
                        normals[3 * id + 2] = cos(phi);
                        texcoords[2 * id] = 1.f - k;
                        texcoords[2 * id + 1] = phi / 2.f / PI;
                }
        }

        for(int i = 0; i < r_size; i++)
        {
                indices[3 * i + 0] = 0;
                indices[3 * i + 1] = circulate_int(i + 2, 1, r_size);
                indices[3 * i + 2] = circulate_int(i + 1, 1, r_size);
        }

        for(int i = 0; i < r_size; i++)
        {
                int id = r_size + 2 * h_size * r_size + i;
                indices[3 * id + 0] = 1 + (h_size + 1) * r_size;
                indices[3 * id + 1] = h_size * r_size + circulate_int(i + 1, 1, r_size);
                indices[3 * id + 2] = h_size * r_size + circulate_int(i + 2, 1, r_size);
        }

        for(int i = 0; i < h_size; i++)
        {
                for(int j = 0; j < r_size; j++)
                {
                        int id = r_size + 2 * i * r_size + 2 * j;
                        
                        int v0 = (i + 1) * r_size + circulate_int(j + 1, 1, r_size);
                        int v1 = i * r_size + circulate_int(j + 1, 1, r_size);
                        int v2 = i * r_size + circulate_int(j + 2, 1, r_size);
                        int v3 = (i + 1) * r_size + circulate_int(j + 2, 1, r_size);

                        indices[3 * id + 0] = v0;
                        indices[3 * id + 1] = v1;
                        indices[3 * id + 2] = v2;
                        indices[3 * id + 3] = v2;
                        indices[3 * id + 4] = v3;
                        indices[3 * id + 5] = v0;
                }
        } 

        mesh.vertex_count = vertices_size;
        mesh.triangle_count = triangles_size;
        mesh.vertices = (float *)malloc(vertices_size * 3 * sizeof(float));
        mesh.normals = (float *)malloc(vertices_size * 3 * sizeof(float));
        mesh.texcoords = (float *)malloc(vertices_size * 2 * sizeof(float));
        mesh.indices = (unsigned short *)malloc(triangles_size * 3 * sizeof(unsigned short));

        memcpy(mesh.vertices, vertices, vertices_size * 3 * sizeof(float));
        memcpy(mesh.normals, normals, vertices_size * 3 * sizeof(float));
        memcpy(mesh.texcoords, texcoords, vertices_size * 2 * sizeof(float));
        memcpy(mesh.indices, indices, triangles_size * 3 * sizeof(unsigned short));

        return mesh;
}

#endif