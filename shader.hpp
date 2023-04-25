#ifndef ENGINE_SHADER
#define ENGINE_SHADER

#include "global.hpp"
#include "material.hpp"
#include "matrix.hpp"
#include "texture.hpp"
#include "triangle.hpp"
#include "vec.hpp"

struct vertex_payload
{
    vec3 vertex_pos[3];
    vec3 vertex_normal[3];
    vec2 vertex_texcoords[3];
};

struct fragment_payload
{
    vec3 fragment_pos;
    vec3 fragment_normal;
    vec2 fragment_texcoords;

    float depth;
};

struct Shader
{
    int triangle_count;
    vertex_payload *ver;

    mat4 model_matrix;
    mat4 view_matrix;
    mat4 projection_matrix;
    mat4 screen_to_world_matrix;

    vec3 camera_pos;
    vec3 camera_dir;

    vec3 light_pos;
    vec3 light_dir;
    vec3 light_radiance;

    texture *albedo_map;
    FBO *shadow_map;
    FBO *G_buffer;

    float kd;
    float ks;
    float ka;
    float metallic;
    float roughness;

    Shader(int _triangle_count) : triangle_count(_triangle_count)
    {
        ver = (vertex_payload *)malloc(triangle_count * sizeof(vertex_payload));
    }
};

bool inside_triangle(const vec2 &v0, const vec2 &v1, const vec2 &v2, const vec2 &xy, float &alpha, float &beta, float &gamma)
{
    return true;
}

void update_fragment_payload(vertex_payload ver, fragment_payload &frag, float alpha, float beta, float gamma)
{
    frag.fragment_pos = alpha * ver.vertex_pos[0] + beta * ver.vertex_pos[1] + gamma * ver.vertex_pos[2];
    frag.fragment_normal = alpha * ver.vertex_normal[0] + beta * ver.vertex_normal[1] + gamma * ver.vertex_normal[2];
    frag.fragment_texcoords = alpha * ver.vertex_texcoords[0] + beta * ver.vertex_texcoords[1] + gamma * ver.vertex_texcoords[2];
}

enum FBO_types
{
    color    = 0,
    depth    = 1,
    normal   = 2,
    position = 3,
    shadow   = 4
};

vec3 pbr_shader(Shader &shader, fragment_payload &frag)
{
    return vec3();
}

#endif