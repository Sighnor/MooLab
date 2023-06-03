#ifndef ENGINE_SHADER
#define ENGINE_SHADER

#include "global.hpp"
#include "material.hpp"
#include "mat.hpp"
#include "texture.hpp"
#include "vec.hpp"
//
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
    Material_Type shader_type;

    int triangle_count;
    vertex_payload *ver;

    mat4 model_matrix;
    mat4 view_matrix;
    mat4 projection_matrix;
    mat4 screen_to_world_matrix;

    vec3 camera_pos;
    vec3 camera_dir;
    float z_near;

    vec3 light_pos;
    vec3 light_dir;
    vec3 light_radiance;

    texture *albedo_map;
    texture *BRDFLut;
    texture *EavgLut;
    FBO *shadow_map;
    FBO *G_buffer;

    float kd;
    float ks;
    float ka;
    float metallic;
    float roughness;

    std::function<vec3(vertex_payload)> vertex_shader;
    std::function<vec3(fragment_payload)> fragment_shader;

    void resize(int _triangle_count);
    fragment_payload get_fragment_payload(
        int ver_id, 
        float alpha, 
        float beta, 
        float gamma,
        float depth);
    vec3 pbr_shader(fragment_payload &frag);

    Shader() 
    {
        triangle_count = 0;
    }
};

void Shader::resize(int _triangle_count)
{
    triangle_count = _triangle_count;
    ver = (vertex_payload *)malloc(triangle_count * sizeof(vertex_payload));

    if (_triangle_count == 0 && triangle_count != 0)
    {
        free(ver);
        ver = NULL;
    }
    else if (_triangle_count > 0 && triangle_count == 0)
    {
        ver = (vertex_payload *)malloc(triangle_count * sizeof(vertex_payload));
        assert(ver != NULL);
    }
    else if (_triangle_count > 0 && triangle_count > 0 && _triangle_count != triangle_count)
    {
        ver = (vertex_payload *)(triangle_count * sizeof(vertex_payload));
        assert(ver != NULL);           
    }
}

fragment_payload Shader::get_fragment_payload(
        int ver_id, 
        float alpha, 
        float beta, 
        float gamma,
        float depth)
{
    fragment_payload frag;

    frag.fragment_pos = vec4_to_vec3(standardize(model_matrix * pos_to_vec4(alpha * ver[ver_id].vertex_pos[0] + beta * ver[ver_id].vertex_pos[1] + gamma * ver[ver_id].vertex_pos[2])));
    frag.fragment_normal = vec4_to_vec3(standardize(model_matrix * vec_to_vec4(alpha * ver[ver_id].vertex_normal[0] + beta * ver[ver_id].vertex_normal[1] + gamma * ver[ver_id].vertex_normal[2])));
    frag.fragment_texcoords = alpha * ver[ver_id].vertex_texcoords[0] + beta * ver[ver_id].vertex_texcoords[1] + gamma * ver[ver_id].vertex_texcoords[2];
    frag.depth = depth;

    return frag;
}

vec3 Shader::pbr_shader(fragment_payload &frag)
{
    vec3 albedo = 0.5f; //powv((*albedo_map)(frag.fragment_texcoords), 2.2f);

    vec3 N = normalize(frag.fragment_normal);
    vec3 V = normalize(camera_pos - frag.fragment_pos);
    float NdotV = std::max(dot(N, V), 0.f);
    
    vec3 F0 = vec3(0.5f); 
    F0 = mix(F0, albedo, metallic);

    vec3 Lo = vec3(0.2f);

    vec3 L = normalize(light_pos - frag.fragment_pos);
    vec3 H = normalize(V + L);
    float NdotL = std::max(dot(N, L), 0.f); 
    vec3 radiance = light_radiance;

    float NDF = DistributionGGX(N, H, roughness);   
    float G   = GeometrySmith(N, V, L, roughness); 
    vec3 F = fresnelSchlick(F0, V, N);
        
    vec3 numerator    = NDF * G * F; 
    float denominator = std::max((4.0 * NdotL * NdotV), 0.001);
    vec3 Fmicro = numerator / denominator;

    vec3 Fms = MultiScatterBRDF(BRDFLut, EavgLut, roughness, frag.fragment_texcoords, NdotL, NdotV);

    vec3 BRDF = Fmicro + Fms;
  
    Lo = Lo + BRDF * radiance * NdotL;
    vec3 color = Lo;
  
    color = color / (color + vec3(1.0));
    color = powv(color, 1.0f / 2.2f); 

    return color;
}

bool inside_2dtriangle(const vec2 &p0, const vec2 &p1, const vec2 &p2, const vec2 &p, float &alpha, float &beta, float &gamma)
{
    alpha = ((p1.y - p2.y) * p.x + (p2.x - p1.x) * p.y + p1.x * p2.y - p2.x * p1.y) / ((p1.y - p2.y) * p0.x + (p2.x - p1.x) * p0.y + p1.x * p2.y - p2.x * p1.y);
    beta = ((p2.y - p0.y) * p.x + (p0.x - p2.x) * p.y + p2.x * p0.y - p0.x * p2.y) / ((p2.y - p0.y) * p1.x + (p0.x - p2.x) * p1.y + p2.x * p0.y - p0.x * p2.y);
    gamma = ((p0.y - p1.y) * p.x + (p1.x - p0.x) * p.y + p0.x * p1.y - p1.x * p0.y) / ((p0.y - p1.y) * p2.x + (p1.x - p0.x) * p2.y + p0.x * p1.y - p1.x * p0.y);

    if(alpha >= 0.f && alpha <= 1.f && beta >= 0.f && beta <= 1.f && gamma >= 0.f && gamma <= 1.f)
    {
        return true;
    }

    return false;
}

bool inside_3dtriangle(const vec3 &v0, const vec3 &v1, const vec3 &v2, const vec3 &v, float &alpha, float &beta, float &gamma)
{
    //v0, v1, v2逆时针
    mat3 tri_matrix = inv_mat(get_coordinate_matrix(cross(v1 - v0, v2 - v0), vec3(0.f, 1.f, 0.f)));

    vec2 p = vec3_to_vec2(tri_matrix * v);
    vec2 p0 = vec3_to_vec2(tri_matrix * v0);
    vec2 p1 = vec3_to_vec2(tri_matrix * v1);
    vec2 p2 = vec3_to_vec2(tri_matrix * v2);

    return inside_2dtriangle(p0, p1, p2, p, alpha, beta, gamma);
}

#endif
