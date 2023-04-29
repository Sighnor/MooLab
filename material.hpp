#ifndef ENGINE_MATERIAL
#define ENGINE_MATERIAL

#include "global.hpp"
#include "texture.hpp"

enum Material_Type
{
    Phong       = 0,
    Bling_Phong = 1,
    PBR         = 2,
};

struct Material
{
    float kd;
    float ks;
    float ka;

    float metallic;
    float roughness;
    texture *tex;
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

#endif