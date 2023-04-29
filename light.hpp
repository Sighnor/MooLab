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

    Point_Light(vec3 pos, vec3 dir, vec3 rad) : light_position(pos), light_direction(dir), light_radiance(rad)
    {
        entity = light_cube(get_coordinate_matrix(dir, vec3(0.f, 1.f, 0.f)), light_position);
    }
};

struct Polygon_Light
{
    std::vector<vec3> light_positions;
    vec3 light_direction;
    vec3 light_radiance;
};


#endif