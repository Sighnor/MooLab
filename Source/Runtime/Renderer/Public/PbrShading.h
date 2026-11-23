#ifndef MOOLAB_PBR_SHADING
#define MOOLAB_PBR_SHADING

// 标准单次散射 metallic-roughness 直接光BRDF的CPU数值核；GLSL对应实现见Shaders/pbr_common.glsl，两者须同步修改。
// 纯数值接口：不依赖GL或存量头文件，只读、无共享写状态，可在光栅并行区内调用。
//
//  F0 = (1-m)*0.04 + m*Albedo；F = F0 + (1-F0)*(1-VdotH)^5
//  D = a2/(PI*b*b)，a2 = r^4，b = (1-NdotH^2) + NdotH^2*a2
//  G = G1(NdotV)*G1(NdotL)，G1(x) = x/(x*(1-k)+k)，k = (r+1)^2/8
//  Specular = D*G*F/max(4*NdotL*NdotV, 0.001)；Diffuse = (1-m)*(1-F)*Albedo/PI
//
// 有限输入的余弦钳至[0,1]，Roughness钳至[0.05,1]，Metallic钳至[0,1]；NdotV或NdotL非正时返回零。
// 单次散射近似：粗糙金属偏暗、介电质掠射处可有能量增益，不保证白炉守恒。

#include "RenderTypes.h"

struct FPbrBrdfInput
{
    // 线性空间albedo
    float Albedo[3] = {};
    float Roughness = 0.5f;
    float Metallic = 0.f;
    float NdotV = 0.f;
    float NdotL = 0.f;
    float NdotH = 0.f;
    float VdotH = 0.f;
};

struct FPbrBrdfOutput
{
    float Diffuse[3] = {};
    float Specular[3] = {};
};

FPbrBrdfOutput EvaluatePbrBrdf(const FPbrBrdfInput& Input);

// 程序化地面的CPU数值核（GPU对应Shaders/ground.frag，两者须同步修改）。约定见RenderTypes.h的FGroundParams。

// 由帧参数与原始平面系数 A*x+B*y+C*z+D=0 准备地面参数：平面按double真实长度归一化，基取与n最不平行的
// 世界轴（并列按X、Y、Z）投影得到U，V = U×n；射线方向由View旋转部分的转置与Proj[0]/Proj[5]得到。
// 平面/投影不可用（非有限、法向为零）或相机不在平面正侧时返回false，调用方应跳过地面。颜色/材质字段保持默认。
bool BuildGroundParams(const FFrameParams& Frame, const double Plane[4], float CellSize, int Width, int Height,
                       FGroundParams& Out);

struct FGroundSample
{
    // 正的view depth（即射线参数t）
    float Depth = 0.f;
    float Position[3] = {};
    // 格子坐标（以格长为单位）与偶格覆盖率；覆盖率已含解析盒滤波，足迹非有限或触及平行奇点时为0.5
    float U = 0.f;
    float V = 0.f;
    float Coverage = 0.f;
    // 线性albedo = Odd + (Even - Odd) * Coverage
    float Albedo[3] = {};
};

// 图像像素列Column、行Row（顶部为0）的地面求交与着色前量。仅当相机在正侧、射线朝向平面且
// ZNear < t <= ClearDepth（t有限）时返回true；平行射线直接拒绝，不夹紧分母。
bool EvaluateGroundPixel(const FGroundParams& Ground, const float CameraPos[3], float ZNear, float ClearDepth,
                         int Column, int Row, FGroundSample& Out);

#endif
