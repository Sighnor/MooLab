#ifndef ENGINE_GLOBAL
#define ENGINE_GLOBAL

#include <iostream>
#include <functional>
#include <math.h>
#include <random>

#define PI 3.14159265358979323846f

static inline float deg_rad(const float& deg) { return deg * M_PI / 180.0; }

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