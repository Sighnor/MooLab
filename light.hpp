#ifndef ENGINE_LIGHT
#define ENGINE_LIGHT

#include "mesh.hpp"
#include "texture.hpp"

struct Ray
{
    vec3 origin;
    vec3 direction;

    Ray(vec3 _origin, vec3 _direction) : origin(_origin), direction(_direction) {}
};

struct Point_Light
{
    vec3 light_position;
    vec3 light_direction;
    vec3 light_radiance;

    Mesh entity;
    FBO *shadow_map;
};

struct Polygon_Light
{
    std::vector<vec3> light_positions;
    vec3 light_direction;
    vec3 light_radiance;
};


#endif