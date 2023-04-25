#ifndef ENGINE_MESH
#define ENGINE_MESH

struct Mesh
{
    int vertex_count;
    int triangle_count;
    float* vertices;
    float* normals;
    float* texcoords;
    int* indices;
};

#endif