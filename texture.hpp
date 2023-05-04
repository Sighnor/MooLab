#ifndef ENGINE_TEXTURE
#define ENGINE_TEXTURE

#include <opencv2\opencv.hpp>
#include <opencv2\core\core.hpp>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include <opencv2\objdetect\objdetect.hpp>
#include <opencv2\imgproc\types_c.h>
#include "array.hpp"
#include <assert.h>
#include "vec.hpp"

struct FBO
{
    int rows;
    int cols;
    vec3* colors;
    float* depths;
    inline vec3& getcolor(int i, int j) const { assert(i >= 0 && i < rows && j >= 0 && j < cols); return colors[i * cols + j]; }
    inline float& getdepth(int i, int j) const { assert(i >= 0 && i < rows && j >= 0 && j < cols); return depths[i * cols + j]; }
    void set(vec3 color, float depth)
    {
        for(int i = 0; i < rows * cols; i++)
        {
            colors[i] = color;
            depths[i] = depth;
        }
    }

    FBO(){}
    FBO(int _rows, int _cols) : rows(_rows), cols(_cols)
    {
        colors = (vec3*)malloc(_rows * _cols * sizeof(vec3));
        depths = (float*)malloc(_rows * _cols * sizeof(float));
    }
};

cv::Mat fbo_to_img(FBO *fbo)
{
    cv::Mat img(fbo->rows, fbo->cols, CV_32FC3, fbo->colors);

    img.convertTo(img, CV_8UC3, 1.0f);
    cv::cvtColor(img, img, cv::COLOR_RGB2BGR);
    
    return img;
}

struct texture
{
    int rows;
    int cols;
    vec3* colors;
    inline vec3& getcolor(int i, int j) const { assert(i >= 0 && i < rows && j >= 0 && j < cols); return colors[i * cols + j]; }
    inline vec3& getcolor(vec2 v) const { assert(v.x >= 0 && v.x < rows && v.y >= 0 && v.y < cols); return colors[float_to_int(v.x) * cols + float_to_int(v.y)]; }

    texture(){}
};

texture img_to_tex(const char* name)
{
    texture tex;
    
    cv::Mat img = cv::imread(name);
    img.convertTo(img, CV_32FC3, 1.0f);

    tex.rows = img.rows;
    tex.cols = img.cols;
    tex.colors = (vec3*)malloc(img.rows * img.cols * sizeof(vec3));
    memcpy(tex.colors, (vec3*)img.data, img.rows * img.cols * sizeof(vec3));

    return tex;
}

#endif