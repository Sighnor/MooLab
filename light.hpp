#ifndef ENGINE_LIGHT
#define ENGINE_LIGHT

#include "material.hpp"
#include "model.hpp"
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

    FBO *shadow_map;

    Point_Light(vec3 pos, vec3 dir, vec3 rad) : light_position(pos), light_direction(dir), light_radiance(rad) {}
};

Model get_light_model(Point_Light &light)
{
    Model model = mesh_material_to_model(light_cube(get_coordinate_matrix(light.light_direction, vec3(0.f, 1.f, 0.f)), light.light_position), Material(Light));

    return model;
}

struct Polygon_Light
{
    std::vector<vec3> light_positions;
    vec3 light_direction;
    vec3 light_radiance;
};


#endif