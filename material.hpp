#ifndef ENGINE_MATERIAL
#define ENGINE_MATERIAL

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

#endif