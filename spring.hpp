#ifndef MOOLAB_SPRING
#define MOOLAB_SPRING

#include "quat.hpp"

vec3 decay_spring_implicit_damping_rot(vec3 v1, vec3 v2, float half_life, float t)
{
    return vec3();
}

vec3 decay_spring_implicit_damping_pos(vec3 v1, vec3 v2, float half_life, float t)
{
    return vec3();
}

float line_function(float a)
{
    return 0.5f * a;
}

#endif