#ifndef ENGINE_CAMERA
#define ENGINE_CAMERA

#include "mat.hpp"
#include "texture.hpp"

struct Camera
{
    float fov;
    float width;
    float height;
    float z_near;
    float z_far;
    vec3 position;
    vec3 direction;
    mat3 orientation;
    FBO *fbo;


    Camera(float _fov, float _width, float _height, float _z_near, float _z_far, vec3 _position, mat3 _orientation) :
    fov(_fov), width(_width), height(_height), z_near(_z_near), z_far(_z_far), position(_position), orientation(_orientation) {}

    void set_pos(vec3 pos);
    void set_view(vec3 view);
    void update_orientation(vec3 look_up);
    mat4 get_view_matrix();
    mat4 get_inverse_view_matrix();
    mat4 get_projection_matrix();
    mat4 get_inverse_projection_matrix();
};

void Camera::set_pos(vec3 pos)
{
    position = pos;
}

void Camera::set_view(vec3 view)
{
    direction = view;
}

void Camera::update_orientation(vec3 look_up)
{
    orientation = get_coordinate_matrix(-direction, look_up);
}

mat4 Camera::get_view_matrix()
{
    return get_inverse_rigid_body_motion_matrix(orientation, position);
}

mat4 Camera::get_inverse_view_matrix()
{
    return get_rigid_body_motion_matrix(orientation, position);
}
//z_far->1, z_near->-1
mat4 Camera::get_projection_matrix()
{
    return get_camera_projection_matrix(z_near, z_far, (width / 2), (height / 2));
}

mat4 Camera::get_inverse_projection_matrix()
{
    return get_inverse_camera_projection_matrix(z_near, z_far, (width / 2), (height / 2));
}

#endif