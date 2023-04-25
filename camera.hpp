#ifndef ENGINE_CAMERA
#define ENGINE_CAMERA

#include "matrix.hpp"
#include "texture.hpp"

struct Camera
{
    int fov;
    int width;
    int height;
    float z_near;
    float z_far;
    vec3 position;
    vec3 direction;
    mat3 orientation;
    FBO *fbo;


    Camera(int _fov, int _width, int _height, float _z_near, float _z_far, vec3 _position, mat3 _orientation) :
    fov(_fov), width(_width), height(_height), z_near(_z_near), z_far(_z_far), position(_position), orientation(_orientation) {}

    void set_view(vec3 view, vec3 look_up);
    mat4 get_view_matrix();
    mat4 get_projection_matrix();
};

void Camera::set_view(vec3 view, vec3 look_up)
{
    direction = view;

    vec3 z_axis = normalize(-view);
    vec3 y_axis = normalize(look_up - dot(look_up, z_axis) * z_axis);
    vec3 x_axis = normalize(cross(y_axis, z_axis));

    orientation.X = x_axis;
    orientation.Y = y_axis;
    orientation.Z = z_axis;
}

mat4 Camera::get_view_matrix()
{
    return mat4(inv_mat(orientation), vec3()) * mat4(eye3(), vec3(-position)) ;
}

mat4 Camera::get_projection_matrix()
{
    return mat4(eye3(), vec3());
}

#endif