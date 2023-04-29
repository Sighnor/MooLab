#ifndef ENGINE_SHADER
#define ENGINE_SHADER

#include "global.hpp"
#include "material.hpp"
#include "mat.hpp"
#include "texture.hpp"
#include "triangle.hpp"
#include "vec.hpp"

enum FBO_types
{
    color    = 0,
    depth    = 1,
    normal   = 2,
    position = 3,
    shadow   = 4
};

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

    vec3 *camera_pos;
    vec3 *camera_dir;

    vec3 *light_pos;
    vec3 *light_dir;
    vec3 *light_radiance;

    texture *albedo_map;
    FBO *shadow_map;
    FBO *G_buffer;

    float kd;
    float ks;
    float ka;
    float metallic;
    float roughness;

    std::function<vec3(vertex_payload)> vertex_shader;
    std::function<vec3(fragment_payload)> fragment_shader;

    vec3 pbr_shader(fragment_payload &frag);

    Shader(int _triangle_count) : triangle_count(_triangle_count)
    {
        ver = (vertex_payload *)malloc(triangle_count * sizeof(vertex_payload));
    }
};

vec3 Shader::pbr_shader(fragment_payload &frag)
{
    vec3 N = normalize(frag.fragment_normal);
    vec3 V = normalize(*camera_pos - frag.fragment_pos);
    float NdotV = std::max(dot(N, V), 0.f);
    
    vec3 F0 = vec3(0.5f); 

    vec3 Lo = vec3(0.1f);

    vec3 L = normalize(*light_pos - frag.fragment_pos);
    vec3 H = normalize(V + L);
    float NdotL = std::max(dot(N, L), 0.f); 

    vec3 radiance = *light_radiance;

    float NDF = DistributionGGX(N, H, roughness);   
    float G   = GeometrySmith(N, V, L, roughness); 
    vec3 F = fresnelSchlick(F0, V, N);
        
    vec3 numerator    = NDF * G * F; 
    float denominator = std::max((4.0 * NdotL * NdotV), 0.001);
    vec3 BRDF = numerator / denominator;

    // print(numerator);

    vec3 color = Lo + BRDF * radiance * NdotL;

    return color;
}

bool inside_2dtriangle(const vec2 &p0, const vec2 &p1, const vec2 &p2, const vec2 &p, float &alpha, float &beta, float &gamma)
{
    alpha = ((p1.y - p2.y) * p.x + (p2.x - p1.x) * p.y + p1.x * p2.y - p2.x * p1.y) / ((p1.y - p2.y) * p0.x + (p2.x - p1.x) * p0.y + p1.x * p2.y - p2.x * p1.y);
    beta = ((p2.y - p0.y) * p.x + (p0.x - p2.x) * p.y + p2.x * p0.y - p0.x * p2.y) / ((p2.y - p0.y) * p1.x + (p0.x - p2.x) * p1.y + p2.x * p0.y - p0.x * p2.y);
    gamma = ((p0.y - p1.y) * p.x + (p1.x - p0.x) * p.y + p0.x * p1.y - p1.x * p0.y) / ((p0.y - p1.y) * p2.x + (p1.x - p0.x) * p2.y + p0.x * p1.y - p1.x * p0.y);

    if(alpha > 0.f && alpha < 1.f && beta > 0.f && beta < 1.f && gamma > 0.f && gamma < 1.f)
    {
        return true;
    }

    return false;
}

bool inside_3dtriangle(const vec3 &v0, const vec3 &v1, const vec3 &v2, const vec3 &v, float &alpha, float &beta, float &gamma)
{
    //v0, v1, v2逆时针
    mat3 tri_matrix = inv_mat(get_coordinate_matrix(cross(v1 - v0, v2 - v0), vec3(0.f, 1.f, 0.f)));

    // print(tri_matrix * inv_mat(tri_matrix));
    // print(tri_matrix);
    // print(inv_mat(tri_matrix));

    vec2 p = vec3_to_vec2(tri_matrix * v);
    vec2 p0 = vec3_to_vec2(tri_matrix * v0);
    vec2 p1 = vec3_to_vec2(tri_matrix * v1);
    vec2 p2 = vec3_to_vec2(tri_matrix * v2);

    return inside_2dtriangle(p0, p1, p2, p, alpha, beta, gamma);
}

void update_fragment_payload(vertex_payload &ver, fragment_payload &frag, float alpha, float beta, float gamma)
{
    frag.fragment_pos = alpha * ver.vertex_pos[0] + beta * ver.vertex_pos[1] + gamma * ver.vertex_pos[2];
    frag.fragment_normal = alpha * ver.vertex_normal[0] + beta * ver.vertex_normal[1] + gamma * ver.vertex_normal[2];
    frag.fragment_texcoords = alpha * ver.vertex_texcoords[0] + beta * ver.vertex_texcoords[1] + gamma * ver.vertex_texcoords[2];
}

#endif