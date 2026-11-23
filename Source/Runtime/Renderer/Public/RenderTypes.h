#ifndef MOOLAB_RENDER_TYPES
#define MOOLAB_RENDER_TYPES

// GPU渲染路径的纯数据与上下文契约：SandboxScene.cpp（存量数学/网格类型）与GL编译单元之间只通过这里交换数据。
// 本头文件不得包含存量头文件（vec.hpp等含非inline定义）、glad或GLFW。
//
// 约定：
//  - 矩阵为16个float列主序（与存量mat4的X/Y/Z/W列布局一致），上传时不转置。
//  - 投影为标准GL形式（near->-1, far->1），不使用glClipControl；深度测试不使用NDC深度，见FFrameParams::ClearDepth。
//  - 图像行0为图像顶部（与OpenCV/imread一致）；GL帧缓冲行0为底部，读回时由渲染器翻转。
//  - 纹理按CPU存储原样上传（行0=图像首行），着色器按CPU的行列寻址规则用texelFetch取值，不做翻转。

#include <cstdint>
#include <vector>

struct FMeshId
{
    int Index = -1;
};

struct FTextureId
{
    int Index = -1;
};

// 网格数据视图（不持有内存），各属性按顶点紧排：Positions/Normals为3个float，TexCoords为2个float
struct FMeshDataView
{
    const float* Positions = nullptr;
    const float* Normals = nullptr;
    const float* TexCoords = nullptr;
    const uint16_t* Indices = nullptr;
    int VertexCount = 0;
    int TriangleCount = 0;
};

// RGB浮点图像视图（不持有内存），取值0~255，与存量texture::colors布局一致
struct FImageRGBFView
{
    const float* Data = nullptr;
    int Width = 0;
    int Height = 0;
};

// 8位RGB图像，行0为顶部，每行Width*3字节紧排
struct FImageRGB8
{
    int Width = 0;
    int Height = 0;
    std::vector<uint8_t> Pixels;
};

struct FFrameParams
{
    float View[16] = {};
    float Proj[16] = {};
    float CameraPos[3] = {};
    float LightPos[3] = {};
    float LightRadiance[3] = {};
    // 清屏色，取值0~1
    float ClearColor[3] = {};
    // 深度沿用CPU参考的定义：观察空间距离（-z）在屏幕空间线性插值；只接受 ZNear < depth <= 已存深度，
    // 深度缓冲清为ClearDepth（CPU为100，即更远的片元一律不绘制）
    float ZNear = 0.f;
    float ClearDepth = 100.f;
};

struct FPbrParams
{
    float Roughness = 0.5f;
    float Metallic = 0.f;
    FTextureId AlbedoMap;
};

// 世界固定的程序化棋盘地面（纯数据，不持有纹理），CPU地面pass与GPU ground.frag按同一规则消费。
// 图像像素列i、行j（顶部为0）的世界射线方向 q = RayTopLeft + i*RayStepX + j*RayStepY，起点为FFrameParams::CameraPos；
// q未归一化且观察空间z分量为-1，故平面求交参数t即正的view depth。采样沿用模型光栅的一行上移：
// ndc = (2*(i+0.5)/W-1, 1-2*(j+1.5)/H)。
struct FGroundParams
{
    // 射线像素增量以此尺寸准备；DrawGround要求与渲染目标一致，禁止静默错位
    int TargetWidth = 0;
    int TargetHeight = 0;
    // 单位法向n与d：n·P + d = 0，相机位于n正侧
    float Plane[4] = {0.f, 1.f, 0.f, 0.f};
    // 格子原点 -d*n 与平面内单位正交基（只依赖平面）
    float Origin[3] = {};
    float AxisU[3] = {1.f, 0.f, 0.f};
    float AxisV[3] = {0.f, 0.f, 1.f};
    float RayTopLeft[3] = {};
    float RayStepX[3] = {};
    float RayStepY[3] = {};
    float CellSize = 0.5f;
    // 线性albedo：(floor(u)+floor(v))为偶数的格取Even
    float AlbedoEven[3] = {0.5f, 0.5f, 0.5f};
    float AlbedoOdd[3] = {0.05f, 0.05f, 0.05f};
    float Roughness = 0.8f;
    float Metallic = 0.f;
};

// 上下文与显示句柄来自Core，不反向依赖窗口。
#include "GlContext.h"
#include "ViewportTexture.h"

enum class EShading : uint8_t
{
    Light,
    Pbr,
};

#endif
