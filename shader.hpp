#ifndef MOOLAB_SHADER
#define MOOLAB_SHADER

#include "global.hpp"
#include "material.hpp"
#include "mat.hpp"
#include "texture.hpp"
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

struct MooShader
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
    vec3 texture_shader(fragment_payload &frag);
    vec3 normal_shader(fragment_payload &frag);

    MooShader() 
    {
        triangle_count = 0;
    }
};

void MooShader::resize(int _triangle_count)
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

fragment_payload MooShader::get_fragment_payload(
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

float DistributionGGX(vec3 N, vec3 H, float roughness)
{
  // TODO: To calculate GGX NDF here
  float a = roughness * roughness;
  float a2 = a * a;
  float NdotH = std::max(dot(N, H), 0.f);
  float NdotH2 = NdotH * NdotH;

  float nom = a2;
  float denom  = (NdotH2 * (a2 - 1.f) + 1.f);
  denom = PI * denom * denom;

  return nom / denom;
}

float GeometrySchlickGGX(float NdotV, float roughness)
{
  // TODO: To calculate Schlick G1 here
  // TODO: To calculate Smith G1 here
  float r = (roughness + 1.f);
  float k = (r * r) / 8.f; 
  float nom = NdotV;
  float denom = NdotV * (1.f - k) + k;

  return nom / denom;
}

float GeometrySmith(vec3 N, vec3 V, vec3 L, float roughness)
{
    // TODO: To calculate Smith G here
    float NdotV = std::max(dot(N, V), 0.f);
    float NdotL = std::max(dot(N, L), 0.f);
    float ggx1 = GeometrySchlickGGX(NdotV, roughness);
    float ggx2 = GeometrySchlickGGX(NdotL, roughness);

    return ggx1 * ggx2;
}

vec3 fresnelSchlick(vec3 F0, vec3 V, vec3 N)
{
    // TODO: To calculate Schlick F here
    float NdotV = dot(N, V);
    return F0 + (vec3(1.f) - F0) * pow(1.f - NdotV, 5.f);
}

vec3 AverageFresnel(vec3 r, vec3 g)
{
    return vec3(0.087237) + 0.0230685*g - 0.0864902*g*g + 0.0774594*g*g*g
           + 0.782654*r - 0.136432*r*r + 0.278708*r*r*r
           + 0.19744*g*r + 0.0360605*g*g*r - 0.2586*g*r*r;
}

vec3 MultiScatterBRDF(texture *albedo_map, texture *BRDFLut, texture *EavgLut, float Roughness, vec2 tex_coords, float NdotL, float NdotV)
{
  vec3 albedo = powv(albedo_map->getcolor(tex_coords) / 255.f, 2.2f);
  // vec3 albedo = albedo_map->getcolor(tex_coords) / 255.f;

  vec3 E_o = BRDFLut->getcolor(NdotL, Roughness) / 255.f;
  vec3 E_i = BRDFLut->getcolor(NdotV, Roughness) / 255.f;

  vec3 E_avg = EavgLut->getcolor(0, Roughness) / 255.f;
  // copper
  vec3 edgetint = vec3(0.827, 0.792, 0.678);
  vec3 F_avg = AverageFresnel(albedo, edgetint);
  
  // TODO: To calculate fms and missing energy here
  vec3 F_add = F_avg * E_avg / (vec3(1.0) - F_avg * (vec3(1.0) - E_avg));
  vec3 F_ms = (vec3(1.0) - E_o) * (vec3(1.0) - E_i) / PI / (vec3(1.0) - E_avg);

  return F_add * F_ms;
}

vec3 MooShader::pbr_shader(fragment_payload &frag)
{
    vec3 albedo = powv(albedo_map->getcolor(frag.fragment_texcoords) / 255.f, 2.2f);
    // vec3 albedo = albedo_map->getcolor(frag.fragment_texcoords) / 255.f;

    vec3 N = normalize(frag.fragment_normal);
    vec3 V = normalize(camera_pos - frag.fragment_pos);
    float NdotV = std::max(dot(N, V), 0.f);
    
    vec3 F0 = vec3(0.5f); 
    F0 = mix(F0, albedo, metallic);

    vec3 Lo = vec3(0.1f);

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

    vec3 Fms = MultiScatterBRDF(albedo_map, BRDFLut, EavgLut, roughness, frag.fragment_texcoords, NdotL, NdotV);

    vec3 BRDF = Fmicro + Fms;
  
    Lo = Lo + BRDF * radiance * NdotL;
    vec3 color = Lo;
  
    color = color / (color + vec3(1.0));
    color = powv(color, 1.0f / 2.2f);

    return color;
}

vec3 MooShader::texture_shader(fragment_payload &frag)
{
    vec3 color = albedo_map->getcolor(frag.fragment_texcoords) / 255.f;
    // vec3 color = vec3(frag.fragment_texcoords.x, frag.fragment_texcoords.y, 0.f);

    return color;
}

vec3 MooShader::normal_shader(fragment_payload &frag)
{
    vec3 color = normalize(frag.fragment_normal);

    return color;
}

#endif