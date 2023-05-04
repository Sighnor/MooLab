#ifndef ENGINE_VECTOR
#define ENGINE_VECTOR

#include "global.hpp"

struct vec2
{
    vec2() : x(0.0f), y(0.0f) {}
    vec2(float a) : x(a), y(a) {}
    vec2(float _x, float _y) : x(_x), y(_y) {}
    
    float x, y;
};

static inline void print(const vec2 &v)
{
    printf("\n%f\n%f\n", v.x, v.y);
}

static inline vec2 operator + (vec2 v1, vec2 v2)
{
    return vec2(v1.x + v2.x, v1.y + v2.y);
}

static inline vec2 operator - (vec2 v1, vec2 v2)
{
    return vec2(v1.x - v2.x, v1.y - v2.y);
}

static inline vec2 operator * (float s, vec2 v)
{
    return vec2(v.x * s, v.y * s);
}

static inline vec2 operator * (vec2 v, float s)
{
    return vec2(v.x * s, v.y * s);
}

static inline vec2 operator / (vec2 v, float s)
{
    return vec2(v.x / s, v.y / s);
}

static inline vec2 operator - (vec2 v)
{
    return vec2(-v.x, -v.y);
}

static inline float dot(vec2 v1, vec2 v2)
{
    return v1.x * v2.x + v1.y * v2.y;
}

static inline float length(const vec2 &v)
{
    return sqrtf(dot(v, v));
}

static inline vec2 normalize(vec2 v, float eps = 1e-8f)
{
    return v / (length(v) + eps);
}

static inline vec2 lerp(vec2 v1, vec2 v2, float alpha)
{
    return v1 * (1.0f - alpha) + v2 * alpha;
}

struct vec3
{
    vec3() : x(0.0f), y(0.0f), z(0.0f) {}
    vec3(float a) : x(a), y(a), z(a) {}
    vec3(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
    
    float x, y, z;
};

static inline void print(const vec3 &v)
{
    printf("\n%f\n%f\n%f\n", v.x, v.y, v.z);
}

static inline vec2 vec3_to_vec2(vec3 v)
{
    return vec2(v.x, v.y);
}

static inline vec3 operator + (vec3 v1, vec3 v2)
{
    return vec3(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

static inline vec3 operator - (vec3 v1, vec3 v2)
{
    return vec3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

static inline vec3 operator * (float s, vec3 v)
{
    return vec3(v.x * s, v.y * s, v.z * s);
}

static inline vec3 operator * (vec3 v, float s)
{
    return vec3(v.x * s, v.y * s, v.z * s);
}

static inline vec3 operator * (vec3 v1, vec3 v2)
{
    return vec3(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
}

static inline vec3 operator / (vec3 v, float s)
{
    return vec3(v.x / s, v.y / s, v.z / s);
}

static inline vec3 operator / (vec3 v1, vec3 v2)
{
    return vec3(v1.x / v2.x, v1.y / v2.y, v1.z / v2.z);
}

static inline vec3 operator - (vec3 v)
{
    return vec3(-v.x, -v.y, -v.z);
}

static inline float dot(vec3 v1, vec3 v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

static inline vec3 cross(vec3 v1, vec3 v2)
{
    return vec3(
        v1.y * v2.z - v1.z * v2.y,
        v1.z * v2.x - v1.x * v2.z,
        v1.x * v2.y - v1.y * v2.x);
}

static inline float length(const vec3 &v)
{
    return sqrtf(dot(v, v));
}

static inline vec3 normalize(vec3 v, float eps = 1e-8f)
{
    return v / (length(v) + eps);
}

static inline vec3 lerp(vec3 v1, vec3 v2, float alpha)
{
    return v1 * (1.0f - alpha) + v2 * alpha;
}

static inline vec3 clampv(vec3 v, float min = 0.f, float max = 255.f)
{
    v.x = clampf(v.x, min, max);
    v.y = clampf(v.y, min, max);
    v.z = clampf(v.z, min, max);

    return v;
}

static inline vec3 powv(vec3 v, float f)
{
    return vec3(pow(v.x, f), pow(v.y, f), pow(v.z, f));
}

static inline vec3 mix(vec3 F0, vec3 albedo, float matallic)
{
    return F0;
}

struct vec4
{
    vec4() : x(), y(), z(), w() {}
    vec4(float a) : x(a), y(a), z(a), w(a) {}
    vec4(float _x, float _y, float _z, float _w) : x(_x), y(_y), z(_z), w(_w) {}
    vec4(vec3 v, float _w) : x(v.x), y(v.y), z(v.z), w(_w) {}
    
    float x, y, z, w;
};

static inline void print(const vec4 &v)
{
    printf("\n%f\n%f\n%f\n%f\n", v.x, v.y, v.z, v.w);
}

static inline vec4 vec_to_vec4(vec3 v)
{
    return vec4(v, 0.0f);
}

static inline vec4 pos_to_vec4(vec3 p)
{
    return vec4(p, 1.0f);
}

static inline vec2 vec4_to_vec2(vec4 v)
{
    return vec2(v.x, v.y);
}

static inline vec3 vec4_to_vec3(vec4 v)
{
    return vec3(v.x, v.y, v.z);
}

static inline vec4 operator + (vec4 v1, vec4 v2)
{
    return vec4(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z, v1.w + v2.w);
}

static inline vec4 operator - (vec4 v1, vec4 v2)
{
    return vec4(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z, v1.w - v2.w);
}

static inline vec4 operator * (float s, vec4 v)
{
    return vec4(v.x * s, v.y * s, v.z * s, v.w * s);
}

static inline vec4 operator * (vec4 v, float s)
{
    return vec4(v.x * s, v.y * s, v.z * s, v.w * s);
}

static inline vec4 operator / (vec4 v, float s)
{
    return vec4(v.x / s, v.y / s, v.z / s, v.w / s);
}

static inline vec4 operator - (vec4 v)
{
    return vec4(-v.x, -v.y, -v.z, -v.w);
}

static inline vec4 normalize(vec4 v, float eps = 1e-8f)
{
    return v / (v.w + eps);
}

#endif