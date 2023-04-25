#ifndef ENGINE_TEXTURE
#define ENGINE_TEXTURE

#include "array.hpp"
#include <assert.h>
#include "vec.hpp"

struct FBO
{
    array2d<vec3> colors;
    inline vec3& operator()(int i, int j) const { assert(i >= 0 && i < colors.rows && j >= 0 && j < colors.cols); return colors(i, j); }
};

struct texture
{
    array2d<vec3> colors;
    inline vec3& operator()(int i, int j) const { assert(i >= 0 && i < colors.rows && j >= 0 && j < colors.cols); return colors(i, j); }
};


#endif