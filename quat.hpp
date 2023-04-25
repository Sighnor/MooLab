#ifndef ENGINE_QUATERNION
#define ENGINE_QUATERNION

#include "vec.hpp"

struct quat
{
    quat() : w(1.f), x(0.f), y(0.f), z(0.f) {}
    quat(float _w, float _x, float _y, float _z) : w(_w), x(_x), y(_y), z(_z) {}
    
    float w, x, y, z;
};

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
    q1.w * q2.x + q1.x * q2.w - q1.y * q2.z + q1.z * q2.y,
    q1.w * q2.y + q1.x * q2.z + q1.y * q2.w - q1.z * q2.x,
    q1.w * q2.z - q1.x * q2.y + q1.y * q2.x + q1.z * q2.w);
}

static inline vec3 operator * (quat q, vec3 v)
{
    return v;
}

static inline quat inv_quat(quat q)
{
    return q;
}

static inline quat rot_vec_to_quat(vec3 rot)
{
    return quat();
}

static inline vec3 quat_to_euler_YXZ(quat q)
{
    return vec3();
}

static inline quat euler_YXZ_to_quat(float x, float y, float z)
{
    return quat();
}

static inline quat slerp(quat q1, quat q2, float alpha)
{
    return quat();
}

static inline vec3 quat_to_avel(quat q)
{
    return vec3();
}

static inline quat avel_to_quat(vec3 v)
{
    return quat();
}


#endif