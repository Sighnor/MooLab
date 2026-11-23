#ifndef MOOLAB_MATERIAL
#define MOOLAB_MATERIAL

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

  //原先除material_type外全部未初始化（metallic直接参与F0混合，结果取决于栈上残值），现给出确定默认值：非金属
  MooMaterial(Material_Type _material_type) : material_type(_material_type), kd(0.f), ks(0.f), ka(0.f), metallic(0.f), roughness(0.5f) {}
};

#endif