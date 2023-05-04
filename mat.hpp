#ifndef ENGINE_MATRIX
#define ENGINE_MATRIX

#include "global.hpp"
#include "vec.hpp"

struct mat3
{
    mat3() : X(), Y(), Z() {}
    mat3(vec3 _X, vec3 _Y, vec3 _Z) : X(_X), Y(_Y), Z(_Z) {}

    vec3 X, Y, Z;
};

static inline void print(const mat3 &m)
{
    printf("\n%f, %f, %f\n", m.X.x, m.Y.x, m.Z.x);
    printf("%f, %f, %f\n", m.X.y, m.Y.y, m.Z.y);
    printf("%f, %f, %f\n", m.X.z, m.Y.z, m.Z.z);
}

static inline mat3 eye3()
{
    return mat3(vec3(1.0f, 0.0f, 0.0f), vec3(0.0f, 1.0f, 0.0f), vec3(0.0f, 0.0f, 1.0f));
}

static inline mat3 operator + (mat3 m1, mat3 m2)
{
    return mat3(m1.X + m2.X, m1.Y + m2.Y, m1.Z + m2.Z);
}

static inline mat3 operator - (mat3 m1, mat3 m2)
{
    return mat3(m1.X - m2.X, m1.Y - m2.Y, m1.Z - m2.Z);
}

static inline mat3 operator * (float s, mat3 m)
{
    return mat3(m.X * s, m.Y * s, m.Z * s);
}

static inline mat3 operator * (mat3 m, float s)
{
    return mat3(m.X * s, m.Y * s, m.Z * s);
}

static inline vec3 operator * (mat3 m, vec3 v)
{
    return vec3(
        (m.X.x * v.x + m.Y.x * v.y + m.Z.x * v.z),
        (m.X.y * v.x + m.Y.y * v.y + m.Z.y * v.z),    
        (m.X.z * v.x + m.Y.z * v.y + m.Z.z * v.z));
}

static inline mat3 operator * (mat3 m1, mat3 m2)
{
    return mat3((m1 * m2.X), (m1 * m2.Y), (m1 * m2.Z));
}

static inline mat3 operator / (mat3 m, float s)
{
    return mat3(m.X / s, m.Y / s, m.Z / s);
}

static inline mat3 vec_to_cross_matrix(vec3 v)
{
    return mat3(
        vec3(0.0f, v.z, -v.y), vec3(-v.z, 0.0f, v.x), vec3(v.y, -v.x, 0.0f));
}

static inline mat3 Rodrigues(float theta, vec3 v)
{
    vec3 axis = normalize(v);
    float phi = deg_to_rad(theta);
    mat3 I = eye3();
    mat3 U = vec_to_cross_matrix(axis);

    return (I + sin(phi) * U + (1.f - cos(phi)) * U * U);
}

static inline float determine(const mat3 &m)
{
    return (m.X.x * (m.Y.y * m.Z.z - m.Y.z * m.Z.y)
          + m.Y.x * (m.Z.y * m.X.z - m.Z.z * m.X.y)
          + m.Z.x * (m.X.y * m.Y.z - m.X.z * m.Y.y));
}

static inline mat3 inv_mat(mat3 m, float eps = 1e-8f)
{
    float det = determine(m);

    if(det < eps)
    {
        return m;
    }

    vec3 X = vec3(
                m.Y.y * m.Z.z - m.Y.z * m.Z.y,
                m.Z.y * m.X.z - m.Z.z * m.X.y,
                m.X.y * m.Y.z - m.X.z * m.Y.y);

    vec3 Y = vec3(
                m.Y.z * m.Z.x - m.Y.x * m.Z.z,
                m.Z.z * m.X.x - m.Z.x * m.X.z,
                m.X.z * m.Y.x - m.X.x * m.Y.z);

    vec3 Z = vec3(
                m.Y.x * m.Z.y - m.Y.y * m.Z.x,
                m.Z.x * m.X.y - m.Z.y * m.X.x,
                m.X.x * m.Y.y - m.X.y * m.Y.x);

    return mat3(X, Y, Z) / det;
}

static inline mat3 get_coordinate_matrix(vec3 z, vec3 y)
{
    vec3 my_z = normalize(z);
    vec3 my_y = normalize(y - dot(y, my_z) * my_z);
    vec3 my_x = normalize(cross(my_y, my_z));

    return mat3(my_x, my_y, my_z);
}

struct mat4
{
    mat4() : X(), Y(), Z(), W() {}
    mat4(float k) : X(k, 0.f, 0.f, 0.f), Y(0.f, k, 0.f, 0.f), Z(0.f, 0.f, k, 0.f), W(0.f, 0.f, 0.f, 1.f) {}
    mat4(vec3 _X, vec3 _Y, vec3 _Z, vec3 _W) : X(vec_to_vec4(_X)), Y(vec_to_vec4(_Y)), Z(vec_to_vec4(_Z)), W(pos_to_vec4(_W)) {}
    mat4(mat3 _Q, vec3 _W) : X(vec_to_vec4(_Q.X)), Y(vec_to_vec4(_Q.Y)), Z(vec_to_vec4(_Q.Z)), W(pos_to_vec4(_W)) {}
    mat4(vec4 _X, vec4 _Y, vec4 _Z, vec4 _W) : X(_X), Y(_Y), Z(_Z), W(_W) {}

    vec4 X, Y, Z, W;
};

static inline void print(const mat4 &m)
{
    printf("\n%f, %f, %f, %f\n", m.X.x, m.Y.x, m.Z.x, m.W.x);
    printf("%f, %f, %f, %f\n", m.X.y, m.Y.y, m.Z.y, m.W.y);
    printf("%f, %f, %f, %f\n", m.X.z, m.Y.z, m.Z.z, m.W.z);
    printf("%f, %f, %f, %f\n", m.X.w, m.Y.w, m.Z.w, m.W.w);
}

static inline mat4 eye4()
{
    return mat4(vec4(1.0f, 0.0f, 0.0f, 0.0f), vec4(0.0f, 1.0f, 0.0f, 0.0f), vec4(0.0f, 0.0f, 1.0f, 0.0f), vec4(0.0f, 0.0f, 0.0f, 1.0f));
}

static inline vec4 operator * (mat4 m, vec4 v)
{
    return vec4(
        (m.X.x * v.x + m.Y.x * v.y + m.Z.x * v.z + m.W.x * v.w),
        (m.X.y * v.x + m.Y.y * v.y + m.Z.y * v.z + m.W.y * v.w),    
        (m.X.z * v.x + m.Y.z * v.y + m.Z.z * v.z + m.W.z * v.w),
        (m.X.w * v.x + m.Y.w * v.y + m.Z.w * v.z + m.W.w * v.w));
}

static inline mat4 operator * (mat4 m1, mat4 m2)
{
    return mat4(m1 * m2.X, m1 * m2.Y, m1 * m2.Z, m1 * m2.W);
}

static inline mat4 get_rigid_body_motion_matrix(mat3 orientation, vec3 position)
{
    return mat4(eye3(), vec3(position)) * mat4(orientation, vec3());
}

static inline mat4 get_inverse_rigid_body_motion_matrix(mat3 orientation, vec3 position)
{
    return mat4(inv_mat(orientation), vec3()) * mat4(eye3(), vec3(-position));
}

static inline mat4 get_camera_projection_matrix(float near, float far, float r, float t)
{
    vec4 x = vec4(near / r, 0.f, 0.f, 0.f);
    vec4 y = vec4(0.f, near / t, 0.f, 0.f);
    vec4 z = vec4(0.f, 0.f, - (far + near) / (far - near), - 1.f);
    vec4 w = vec4(0.f, 0.f, - 2 * near * far / (far - near), 0.f);
    return mat4(x, y, z, w);
}

static inline mat4 get_inverse_camera_projection_matrix(float near, float far, float r, float t)
{
    vec4 x = vec4(r / near, 0.f, 0.f, 0.f);
    vec4 y = vec4(0.f, t / near, 0.f, 0.f);
    vec4 z = vec4(0.f, 0.f, 0.f, - (far - near) / (2 * near * far));
    vec4 w = vec4(0.f, 0.f, - 1.f, (far + near) / (2 * near * far));
    return mat4(x, y, z, w);
}

#endif