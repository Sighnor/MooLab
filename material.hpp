#ifndef ENGINE_MATERIAL
#define ENGINE_MATERIAL

#include "global.hpp"
#include "texture.hpp"

enum Material_Type
{
  LIGHT       = 0,
  PBR         = 1,
  TEXTURE     = 2,
  NORMALS     = 3,
};

struct MooMaterial
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

  MooMaterial(Material_Type _material_type) : material_type(_material_type) {}
};

#endif