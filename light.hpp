#ifndef ENGINE_LIGHT
#define ENGINE_LIGHT

#include "material.hpp"
#include "model.hpp"
#include "texture.hpp"

struct MooRay
{
    vec3 origin;
    vec3 direction;

    MooRay(vec3 _origin, vec3 _direction) : origin(_origin), direction(_direction) {}
};

struct Point_Light
{
    vec3 light_position;
    vec3 light_direction;
    vec3 light_radiance;

    FBO *shadow_map;

    Point_Light(vec3 pos, vec3 dir, vec3 rad) : light_position(pos), light_direction(dir), light_radiance(rad) {}
};

MooModel get_light_model(Point_Light &light)
{
    MooModel model = mesh_material_to_model(cube(), MooMaterial(LIGHT));
    model.transform = mat4(eye3(), light.light_position) * mat4(0.01f);

    return model;
}

struct Polygon_Light
{
    array1d<vec3> light_positions;
    vec3 light_direction;
    vec3 light_radiance;
};


#endif