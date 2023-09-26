#ifndef MOOLAB_TEXTURE
#define MOOLAB_TEXTURE

#include <opencv2\opencv.hpp>
#include <opencv2\core\core.hpp>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include <opencv2\objdetect\objdetect.hpp>
#include <opencv2\imgproc\types_c.h>
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

FBO concatenate_two_fbos(FBO *fbo0, FBO *fbo1, float r = 1.f, float c = 1.f)
{
    FBO fbo(r * fbo0->rows, c * (fbo0->cols + fbo1->cols));

    for(int i = 0; i < fbo.rows; i++)
    {
        for(int j = 0; j < 0.5f * fbo.cols; j++)
        {
            fbo.getcolor(i, j) = fbo0->getcolor(0.5f * (1.f - r) * fbo0->rows + i, 0.5f * (1.f - c) * fbo0->cols + j);
            fbo.getcolor(i, c * fbo0->cols + j) = fbo1->getcolor(0.5f * (1.f - r) * fbo0->rows + i, 0.5f * (1.f - c) * fbo0->cols + j);

            fbo.getdepth(i, j) = fbo0->getdepth(0.5f * (1.f - r) * fbo0->rows + i, 0.5f * (1.f - c) * fbo0->cols + j);
            fbo.getdepth(i, c * fbo0->cols + j) = fbo1->getdepth(0.5f * (1.f - r) * fbo0->rows + i, 0.5f * (1.f - c) * fbo0->cols + j);
        }
    }

    return fbo;
}

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
    inline vec3& getcolor(float i, float j) const 
    { 
        int y = float_to_int(clampf(i * rows, 0.f, rows - 1));
        int x = float_to_int(clampf(j * cols, 0.f, cols - 1));
        return colors[y * cols + x]; 
    }
    inline vec3& getcolor(vec2 v) const 
    {
        int y = float_to_int(clampf(v.x * rows, 0.f, rows - 1));
        int x = float_to_int(clampf(v.y * cols, 0.f, cols - 1));
        return colors[y * cols + x]; 
    }

    texture(){}
};

texture img_to_tex(const char* name)
{
    texture tex;
    
    cv::Mat img = cv::imread(name);
    cv::cvtColor(img, img, cv::COLOR_BGR2RGB);
    img.convertTo(img, CV_32FC3, 1.0f);

    tex.rows = img.rows;
    tex.cols = img.cols;
    tex.colors = (vec3*)malloc(img.rows * img.cols * sizeof(vec3));
    memcpy(tex.colors, (vec3*)img.data, img.rows * img.cols * sizeof(vec3));

    return tex;
}

cv::Mat tex_to_img(texture *tex)
{
    cv::Mat img(tex->rows, tex->cols, CV_32FC3, tex->colors);

    img.convertTo(img, CV_8UC3, 1.0f);
    cv::cvtColor(img, img, cv::COLOR_RGB2BGR);
    
    return img;
}

#endif