#ifndef ENGINE_MATERIAL
#define ENGINE_MATERIAL

#include "global.hpp"
#include "texture.hpp"

enum Material_Type
{
    Light       = 0,
    PBR         = 1,
};

struct Material
{
  Material_Type material_type;

  float kd;
  float ks;
  float ka;

  float metallic;
  float roughness;
  texture tex;
  texture BRDFLut;
  texture EavgLut;

  Material(Material_Type _material_type) : material_type(_material_type) {}
};

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

vec3 MultiScatterBRDF(texture *BRDFLut, texture *EavgLut, float Roughness, vec2 tex_coords, float NdotL, float NdotV)
{
  vec3 albedo = 0.5f; //powv((*BRDFLut)(tex_coords), 2.2f) / 255.f;

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

#endif