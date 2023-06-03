#ifndef ENGINE_QUATERNION
#define ENGINE_QUATERNION

#include "mat.hpp"

struct quat
{
    quat() : w(1.f), x(0.f), y(0.f), z(0.f) {}
    quat(float _w, float _x, float _y, float _z) : w(_w), x(_x), y(_y), z(_z) {}
    quat(float theta, vec3 axis) 
    {
        axis = normalize(axis);
        float phi = deg_to_rad(theta);
        float cos_phi_2 = cos(phi / 2.f);
        float sin_phi_2 = sin(phi / 2.f);

        w = cos_phi_2;
        x = sin_phi_2 * axis.x;
        y = sin_phi_2 * axis.y;
        z = sin_phi_2 * axis.z;
    }
    
    float w, x, y, z;
};

static inline void print(const quat &q)
{
    printf("\n%f\n%f\n%f\n%f\n", q.w, q.x, q.y, q.z);
}

static inline quat operator + (quat q1, quat q2)
{
    return quat(q1.w + q2.w, q1.x + q2.x, q1.y + q2.y, q1.z + q2.z);
}

static inline quat operator - (quat q)
{
    return quat(-q.w, -q.x, -q.y, -q.z);
}

static inline quat operator * (quat q, float s)
{
    return quat(q.w * s, q.x * s, q.y * s, q.z * s);
}

static inline quat operator * (float s, quat q)
{
    return quat(q.w * s, q.x * s, q.y * s, q.z * s);
}

static inline quat operator * (quat q1, quat q2)
{
    return quat(
    q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z,
    q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y,
    q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x,
    q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w);
}

// quat = cos(t / 2) + sin(t / 2) * R;
// inv_quat = cos(t / 2) - sin(t / 2) * R;
// q * v = v + 2 * cos(t / 2) * sin(t / 2) * cross(R, v) + 2 * sin(t / 2) * sin(t / 2) * cross(R, cross(R, v));
//       = v + sin(t) * cross(R, v) + (1 - cos(t)) * cross(R, cross(R, v));
// - cos(t / 2) + sin(t / 2) * R equal cos(- t / 2) + sin(- t / 2) * R

static inline vec3 operator * (quat q, vec3 v)
{
    vec3 X = vec3(
        1.f - 2.f * q.y * q.y - 2.f * q.z * q.z, 
        2.f * q.w * q.z + 2.f * q.x * q.y, 
        - 2.f * q.w * q.y + 2.f * q.x * q.z);
    vec3 Y = vec3(
        - 2.f * q.w * q.z + 2.f * q.x * q.y,
        1.f - 2.f * q.x * q.x - 2.f * q.z * q.z,
        2.f * q.w * q.x + 2.f * q.y * q.z);
    vec3 Z = vec3(
        2.f * q.w * q.y + 2.f * q.x * q.z,
        - 2.f * q.w * q.x + 2.f * q.y * q.z,
        1.f - 2.f * q.x * q.x - 2.f * q.y * q.y);
    return mat3(X, Y, Z) * v;
}

// static inline vec3 operator * (quat q, vec3 v)
// {
//     quat qv = (
//         quat(q.w, q.x, q.y, q.z) * 
//         quat(0.f, v.x, v.y, v.z) * 
//         quat(q.w, -q.x, -q.y, -q.z));

//     return vec3(qv.x, qv.y, qv.z);
// }

static inline float dot(quat q1, quat q2)
{
    return (q1.w * q2.w + q1.x * q2.x + q1.y * q2.y + q1.z * q2.z);
}

static inline quat inv_quat(quat q)
{
    return quat(q.w, -q.x, -q.y, -q.z);
}

static inline quat rot_vec_to_quat(vec3 rot)
{
    float phi = length(rot);
    if(phi < 0.00000001f)
    {
        return quat(1.f, 0.f, 0.f, 0.f);
    }

    return quat(rad_to_deg(phi), rot);
}

// (1 - 2 * y * y - 2 * z * z       , - 2 * w * z + 2 * x * y         , 2 * w * y + 2 * x * z
//  2 * w * z + 2 * x * y           , 1 - 2 * x * x - 2 * z * z       , - 2 * w * x + 2 * y * z
//  - 2 * w * y + 2 * x * z         , 2 * w * x + 2 * y * z           , 1 - 2 * x * x - 2 * y * y)

// (cosy * cosz + sinx * siny * sinz, cosz * sinx * siny - cosy * sinz, cosx * siny
//  cosx * sinz                     , cosx * cosz                     , - sinx
//  cosy * sinx * sinz - cosz * siny, siny * sinz + cosy * cosz * sinx, cosx * cosy)

//asin - pi / 2 to pi / 2
//acos 0 to pi

static inline vec3 quat_to_euler_YXZ(quat q)
{
    float x, y, z;
    x = asin(2 * q.w * q.x - 2 * q.y * q.z);
    y = asin((2 * q.w * q.y + 2 * q.x * q.z) / cos(x));
    z = asin((2 * q.w * q.z + 2 * q.x * q.y) / cos(x));
    return vec3(rad_to_deg(y), rad_to_deg(x), rad_to_deg(z));
}

static inline quat euler_YXZ_to_quat(float y, float x, float z)
{
    return (
        quat(y, vec3(0.f, 1.f, 0.f)) * 
        quat(x, vec3(1.f, 0.f, 0.f)) * 
        quat(z, vec3(0.f, 0.f, 1.f)));
}

// length(v) = 1
// v1 = q1 * v
// v2 = q2 * v
// phi = arccos(v1 * v2)
// phi1 = (1 - alpha) * phi
// phi2 = alpha * phi
// v1 = cos(phi2) * v3 + sin(phi2) * v4
// v2 = cos(phi1) * v3 - sin(phi1) * v4
// v3 = sin(phi1) / sin(phi1 + phi2) * v1 + sin(phi2) / sin(phi1 + phi2) * v2
//    = sin(phi1) / sin(phi) * v1 + sin(phi2) / sin(phi) * v2
// q3 * v = (sin(phi1) / sin(phi) * q1 + sin(phi2) / sin(phi) * q2) * v
// q3 = sin(phi1) / sin(phi) * q1 + sin(phi2) / sin(phi) * q2

static inline quat slerp(quat q1, quat q2, float alpha)
{
    float scale1;
    float scale2;
    float dot_q1q2 = dot(q1, q2);

    float phi = acos(abs(dot_q1q2));

    if(phi < 1e-4f)
    {
        scale1 = 1.f - alpha;
        scale2 = alpha;
    }
    else
    {
        scale1 = sin((1.f - alpha) * phi) / sin(phi);
        scale2 = sin(alpha * phi) / sin(phi);
    }

    if(dot_q1q2 < 0.f)
    {
        scale2 = -scale2;
    }

    return scale1 * q1 + scale2 * q2;
}

static inline vec3 quat_to_avel(quat last, quat curr, float dt = 0.0166667f)
{
    quat diff = curr * inv_quat(last);
    float scale;
    if(diff.w > 0.99f)
    {
        scale = 2.f / dt;
    }
    else
    {
        float theta = 2.f * acos(diff.w);
        scale = theta / sin(theta / 2.f) / dt;
    }
    return scale * vec3(diff.x, diff.y, diff.z);
}

static inline quat avel_to_quat(vec3 avel, float t)
{
    float theta = rad_to_deg(length(avel) * t);
    return quat(theta, avel);
}


#endif