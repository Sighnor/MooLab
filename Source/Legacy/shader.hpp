#ifndef MOOLAB_SHADER
#define MOOLAB_SHADER

#include "Core.h"
#include "PbrShading.h"
#include "global.hpp"
#include "material.hpp"
#include "mat.hpp"
#include "texture.hpp"
#include "vec.hpp"

#include <type_traits>

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
    Material_Type shader_type = LIGHT;

    int triangle_count = 0;
    vertex_payload *ver = NULL;

    mat4 model_matrix;
    mat4 view_matrix;
    mat4 projection_matrix;
    mat4 screen_to_world_matrix;

    vec3 camera_pos;
    vec3 camera_dir;
    float z_near = 0.f;

    vec3 light_pos;
    vec3 light_dir;
    vec3 light_radiance;

    texture *albedo_map = NULL;
    FBO *shadow_map = NULL;
    FBO *G_buffer = NULL;

    float kd = 0.f;
    float ks = 0.f;
    float ka = 0.f;
    float metallic = 0.f;
    float roughness = 0.f;

    //只允许所有者在按位拷贝出副本之前调用：副本共享ver，任一副本resize都会使其余副本悬空
    void resize(int _triangle_count);
    //以下着色函数在光栅并行区内被多线程调用，必须保持const（不写成员）
    fragment_payload get_fragment_payload(
        int ver_id,
        float alpha,
        float beta,
        float gamma,
        float depth) const;

    vec3 pbr_shader(const fragment_payload &frag) const;
    vec3 texture_shader(const fragment_payload &frag) const;
    vec3 normal_shader(const fragment_payload &frag) const;
};

//MooModel用malloc分配MooShader数组并按位拷贝，故本结构体须保持可平凡拷贝（不得含std::function/std::vector等成员）
static_assert(std::is_trivially_copyable<MooShader>::value, "MooShader must stay trivially copyable");

void MooShader::resize(int _triangle_count)
{
    MOO_CHECK(_triangle_count >= 0);
    //仅在尺寸变化时重新分配；原实现先覆盖triangle_count再判断，导致每次无条件malloc且旧缓冲泄漏
    if(_triangle_count == triangle_count)
    {
        return;
    }

    free(ver);
    ver = NULL;
    triangle_count = std::max(_triangle_count, 0);

    if(triangle_count > 0)
    {
        ver = (vertex_payload *)malloc(triangle_count * sizeof(vertex_payload));
        MOO_CHECK_ALWAYS(ver != NULL, "out of memory allocating vertex payloads");
    }
}

fragment_payload MooShader::get_fragment_payload(
        int ver_id,
        float alpha,
        float beta,
        float gamma,
        float depth) const
{
    fragment_payload frag;

    frag.fragment_pos = vec4_to_vec3(standardize(model_matrix * pos_to_vec4(alpha * ver[ver_id].vertex_pos[0] + beta * ver[ver_id].vertex_pos[1] + gamma * ver[ver_id].vertex_pos[2])));
    frag.fragment_normal = vec4_to_vec3(standardize(model_matrix * vec_to_vec4(alpha * ver[ver_id].vertex_normal[0] + beta * ver[ver_id].vertex_normal[1] + gamma * ver[ver_id].vertex_normal[2])));
    frag.fragment_texcoords = alpha * ver[ver_id].vertex_texcoords[0] + beta * ver[ver_id].vertex_texcoords[1] + gamma * ver[ver_id].vertex_texcoords[2];
    frag.depth = depth;

    return frag;
}

//线性albedo的点光着色（无距离衰减）：常量底光0.1 + BRDF*radiance*NdotL，再Reinhard与gamma 1/2.2，返回[0, 1]颜色。
//网格与程序化地面共用；BRDF见PbrShading.h，GLSL对应MooShadePbr（pbr_common.glsl）。N须已归一化。
vec3 MooShadePbr(const vec3& Albedo, const vec3& Normal, const vec3& FragPos, const vec3& CameraPos,
                 const vec3& LightPos, const vec3& LightRadiance, float Roughness, float Metallic)
{
    const vec3 V = normalize(CameraPos - FragPos);
    const vec3 L = normalize(LightPos - FragPos);
    const vec3 H = normalize(V + L);

    FPbrBrdfInput Input;
    Input.Albedo[0] = Albedo.x;
    Input.Albedo[1] = Albedo.y;
    Input.Albedo[2] = Albedo.z;
    Input.Roughness = Roughness;
    Input.Metallic = Metallic;
    Input.NdotV = dot(Normal, V);
    Input.NdotL = dot(Normal, L);
    Input.NdotH = dot(Normal, H);
    Input.VdotH = dot(V, H);
    const FPbrBrdfOutput Brdf = EvaluatePbrBrdf(Input);
    const vec3 TotalBrdf(Brdf.Diffuse[0] + Brdf.Specular[0], Brdf.Diffuse[1] + Brdf.Specular[1], Brdf.Diffuse[2] + Brdf.Specular[2]);

    const float NdotL = std::max(dot(Normal, L), 0.f);
    vec3 Color = vec3(0.1f) + TotalBrdf * LightRadiance * NdotL;

    Color = Color / (Color + vec3(1.0));
    Color = powv(Color, 1.0f / 2.2f);

    return Color;
}

vec3 MooShader::pbr_shader(const fragment_payload &frag) const
{
    vec3 albedo = powv(albedo_map->getcolor(frag.fragment_texcoords) / 255.f, 2.2f);
    vec3 N = normalize(frag.fragment_normal);

    return MooShadePbr(albedo, N, frag.fragment_pos, camera_pos, light_pos, light_radiance, roughness, metallic);
}

vec3 MooShader::texture_shader(const fragment_payload &frag) const
{
    vec3 color = albedo_map->getcolor(frag.fragment_texcoords) / 255.f;
    // vec3 color = vec3(frag.fragment_texcoords.x, frag.fragment_texcoords.y, 0.f);

    return color;
}

vec3 MooShader::normal_shader(const fragment_payload &frag) const
{
    vec3 color = normalize(frag.fragment_normal);

    return color;
}

#endif
