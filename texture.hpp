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
    array2d<vec3> colors;
    inline vec3& operator()(int i, int j) const { assert(i >= 0 && i < colors.rows && j >= 0 && j < colors.cols); return colors(i, j); }

    FBO(){}
    FBO(int _rows, int _cols, vec3 color) : rows(_rows), cols(_cols)
    {
        colors.resize(rows, cols);
        colors.set(color);
    }
};

struct texture
{
    array2d<vec3> colors;
    inline vec3& operator()(int i, int j) const { assert(i >= 0 && i < colors.rows && j >= 0 && j < colors.cols); return colors(i, j); }
};

cv::Mat FBO_to_IMG(FBO *fbo)
{
    cv::Mat img(fbo->rows, fbo->cols, CV_32FC3, fbo->colors.data);

    img.convertTo(img, CV_8UC3, 1.0f);
    cv::cvtColor(img, img, cv::COLOR_RGB2BGR);
    
    return img;
}


#endif