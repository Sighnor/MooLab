#ifndef ENGINE_GLOBAL
#define ENGINE_GLOBAL

#include <algorithm>
#include <iostream>
#include <functional>
#include <math.h>
#include <map>
#include <random>
#include <vector>

#define PI 3.14159265358979323846f
//弧度、角度转化
static inline float deg_to_rad(const float& deg) { return deg * PI / 180.f; }
static inline float rad_to_deg(const float& rad) { return rad / PI * 180.f; }

static inline float clampf(float x, float min, float max)
{
    if(x < min)
    {
        x = min;
    }
    else if(x > max)
    {
        x = max;
    }
    return x;
}

static inline int float_to_int(float f)
{
    float a = floor(f);
    float b = ceil(f);

    if(f - a > 0.5f)
    {
        return b;
    }
    return a;
}

static inline float maxf_3(float x, float y, float z)
{
    if(x > y)
    {
        if(x > z)
        {
            return x;
        }
        else
        {
            return z;
        }
    }
    else
    {
        if(y > z)
        {
            return y;
        }
        else
        {
            return z;
        }
    }
}

static inline float minf_3(float x, float y, float z)
{
    if(x < y)
    {
        if(x < z)
        {
            return x;
        }
        else
        {
            return z;
        }
    }
    else
    {
        if(y < z)
        {
            return y;
        }
        else
        {
            return z;
        }
    }
}

bool compare_map_int_float(const std::pair<int, float> x, const std::pair<int, float> y)
{
    return x.second < y.second;
}

static inline float circulate_int(int f, int min, int max)
{
    if(f < min)
    {
        return f + (max - min);
    }
    else if(f > max)
    {
        return f - (max - min);
    }
    return f;
}

static inline float circulate_float(float f, float min, float max)
{
    if(f < min)
    {
        return f + (max - min);
    }
    else if(f > max)
    {
        return f - (max - min);
    }
    return f;
}

static inline float get_random_float()
{
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<float> dist(0.0f, 1.f);

    return dist(rng);
}

static inline float get_random_int()
{
    std::random_device rd;
    std::mt19937 eng(rd());
    std::uniform_int_distribution<int> dis(10, 20);

    return dis(eng);
}

#endif
