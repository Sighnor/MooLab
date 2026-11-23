#include "PbrShading.h"

#include <cmath>

namespace
{
constexpr float Pi = 3.14159265358979323846f;
// 足迹宽度下限（格）防float相消；达到上限时该方向取平均色
constexpr float GroundMinFootprint = 1e-3f;
constexpr float GroundMaxFootprint = 64.f;

// 用比较实现钳制，运算顺序与GLSL版逐项一致
float ClampUnit(float Value, float Min, float Max)
{
    if(Value < Min) Value = Min;
    if(Value > Max) Value = Max;
    return Value;
}

float SchlickG1(float X, float K)
{
    return X / (X * (1.f - K) + K);
}

float Dot3(const float A[3], const float B[3])
{
    return A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
}

// 周期2的平方波（floor(x)偶为+1、奇为-1）的原函数
float SquareWaveIntegral(float X)
{
    return 1.f - std::fabs(X - 2.f * std::floor(X / 2.f) - 1.f);
}

// 平方波在[X-W/2, X+W/2]上的盒滤波均值
float FilteredSquareWave(float X, float W)
{
    if(W >= GroundMaxFootprint) return 0.f;
    if(W < GroundMinFootprint) W = GroundMinFootprint;
    return ClampUnit((SquareWaveIntegral(X + 0.5f * W) - SquareWaveIntegral(X - 0.5f * W)) / W, -1.f, 1.f);
}
}

FPbrBrdfOutput EvaluatePbrBrdf(const FPbrBrdfInput& Input)
{
    FPbrBrdfOutput Output;
    // 取反比较使NaN余弦也走零分支
    if(!(Input.NdotV > 0.f) || !(Input.NdotL > 0.f)) return Output;

    const float NdotV = ClampUnit(Input.NdotV, 0.f, 1.f);
    const float NdotL = ClampUnit(Input.NdotL, 0.f, 1.f);
    const float NdotH = ClampUnit(Input.NdotH, 0.f, 1.f);
    const float VdotH = ClampUnit(Input.VdotH, 0.f, 1.f);
    const float Roughness = ClampUnit(Input.Roughness, 0.05f, 1.f);
    const float Metallic = ClampUnit(Input.Metallic, 0.f, 1.f);

    const float Alpha = Roughness * Roughness;
    const float Alpha2 = Alpha * Alpha;
    const float NdotH2 = NdotH * NdotH;
    const float B = (1.f - NdotH2) + NdotH2 * Alpha2;
    const float D = Alpha2 / (Pi * B * B);

    const float RoughnessPlusOne = Roughness + 1.f;
    const float K = RoughnessPlusOne * RoughnessPlusOne / 8.f;
    const float G = SchlickG1(NdotV, K) * SchlickG1(NdotL, K);

    const float X = 1.f - VdotH;
    const float X2 = X * X;
    const float X5 = X2 * X2 * X;
    const float Denominator = 4.f * NdotL * NdotV > 0.001f ? 4.f * NdotL * NdotV : 0.001f;

    for(int Channel = 0; Channel < 3; ++Channel)
    {
        const float Albedo = Input.Albedo[Channel];
        const float F0 = (1.f - Metallic) * 0.04f + Metallic * Albedo;
        const float F = F0 + (1.f - F0) * X5;
        Output.Specular[Channel] = D * G * F / Denominator;
        Output.Diffuse[Channel] = (1.f - Metallic) * (1.f - F) * Albedo / Pi;
    }
    return Output;
}

bool BuildGroundParams(const FFrameParams& Frame, const double Plane[4], float CellSize, int Width, int Height,
                       FGroundParams& Out)
{
    if(Width <= 0 || Height <= 0 || !std::isfinite(CellSize) || !(CellSize > 0.f)) return false;
    const double Length = std::sqrt(Plane[0] * Plane[0] + Plane[1] * Plane[1] + Plane[2] * Plane[2]);
    if(!std::isfinite(Length) || !(Length > 0.) || !std::isfinite(Plane[3])) return false;
    const double Normal[3] = {Plane[0] / Length, Plane[1] / Length, Plane[2] / Length};
    const double Offset = Plane[3] / Length;

    // 与n最不平行的世界轴：|n_i|最小，并列取靠前者
    int Axis = 0;
    for(int Index = 1; Index < 3; ++Index)
        if(std::fabs(Normal[Index]) < std::fabs(Normal[Axis])) Axis = Index;
    double U[3] = {0., 0., 0.};
    U[Axis] = 1.;
    const double Along = U[Axis] * Normal[Axis];
    for(int Index = 0; Index < 3; ++Index) U[Index] -= Along * Normal[Index];
    const double ULength = std::sqrt(U[0] * U[0] + U[1] * U[1] + U[2] * U[2]);
    for(double& Value : U) Value /= ULength;
    const double V[3] = {U[1] * Normal[2] - U[2] * Normal[1], U[2] * Normal[0] - U[0] * Normal[2], U[0] * Normal[1] - U[1] * Normal[0]};

    // 观察空间射线 (ndc.x/Proj[0], ndc.y/Proj[5], -1)，与模型光栅同样采样行中心并整体上移一行
    const double ProjX = Frame.Proj[0];
    const double ProjY = Frame.Proj[5];
    if(!std::isfinite(ProjX) || !std::isfinite(ProjY) || ProjX == 0. || ProjY == 0.) return false;
    const double TopLeftView[3] = {(1. / Width - 1.) / ProjX, (1. - 3. / Height) / ProjY, -1.};
    const double StepXView[3] = {2. / (Width * ProjX), 0., 0.};
    const double StepYView[3] = {0., -2. / (Height * ProjY), 0.};
    // 世界方向 = R^T * v，R为View列主序左上3x3：world_i = sum_k View[4*i+k]*v_k
    auto ToWorld = [&Frame](const double View[3], float Result[3])
    {
        for(int Index = 0; Index < 3; ++Index)
            Result[Index] = static_cast<float>(Frame.View[4 * Index] * View[0] + Frame.View[4 * Index + 1] * View[1]
                                               + static_cast<double>(Frame.View[4 * Index + 2]) * View[2]);
    };

    FGroundParams Result = Out;
    Result.TargetWidth = Width;
    Result.TargetHeight = Height;
    for(int Index = 0; Index < 3; ++Index)
    {
        Result.Plane[Index] = static_cast<float>(Normal[Index]);
        Result.Origin[Index] = static_cast<float>(-Offset * Normal[Index]);
        Result.AxisU[Index] = static_cast<float>(U[Index]);
        Result.AxisV[Index] = static_cast<float>(V[Index]);
    }
    Result.Plane[3] = static_cast<float>(Offset);
    ToWorld(TopLeftView, Result.RayTopLeft);
    ToWorld(StepXView, Result.RayStepX);
    ToWorld(StepYView, Result.RayStepY);
    Result.CellSize = CellSize;

    // 转float后非有限的平面或射线不可渲染，安全跳过
    const float* Fields[] = {Result.Plane, Result.Origin, Result.AxisU, Result.AxisV, Result.RayTopLeft, Result.RayStepX, Result.RayStepY};
    for(const float* Field : Fields)
        for(int Index = 0; Index < 3; ++Index)
            if(!std::isfinite(Field[Index])) return false;
    if(!std::isfinite(Result.Plane[3])) return false;
    const float Side = Dot3(Result.Plane, Frame.CameraPos) + Result.Plane[3];
    if(!(Side > 0.f)) return false;
    Out = Result;
    return true;
}

bool EvaluateGroundPixel(const FGroundParams& Ground, const float CameraPos[3], float ZNear, float ClearDepth,
                         int Column, int Row, FGroundSample& Out)
{
    const float* N = Ground.Plane;
    const float Side = Dot3(N, CameraPos) + Ground.Plane[3];
    if(!(Side > 0.f)) return false;

    const float I = static_cast<float>(Column);
    const float J = static_cast<float>(Row);
    float Ray[3];
    for(int Index = 0; Index < 3; ++Index)
        Ray[Index] = Ground.RayTopLeft[Index] + I * Ground.RayStepX[Index] + J * Ground.RayStepY[Index];
    // 只接受朝向平面的射线；平行（含分母为0）与NaN一并拒绝
    const float NdotRay = Dot3(N, Ray);
    if(!(NdotRay < 0.f)) return false;
    const float T = -Side / NdotRay;
    if(!(T > ZNear) || !(T <= ClearDepth)) return false;

    float Relative[3];
    for(int Index = 0; Index < 3; ++Index)
    {
        Relative[Index] = (CameraPos[Index] - Ground.Origin[Index]) + T * Ray[Index];
        Out.Position[Index] = CameraPos[Index] + T * Ray[Index];
    }
    const float InvCell = 1.f / Ground.CellSize;
    const float U = Dot3(Relative, Ground.AxisU) * InvCell;
    const float V = Dot3(Relative, Ground.AxisV) * InvCell;

    // 交点对像素轴的解析导数：dt = -t*(n·dq)/(n·q)，dP = t*dq + dt*q
    float Footprint[2] = {0.f, 0.f};
    const float* Steps[2] = {Ground.RayStepX, Ground.RayStepY};
    for(const float* Step : Steps)
    {
        const float DeltaT = -T * Dot3(N, Step) / NdotRay;
        float DeltaP[3];
        for(int Index = 0; Index < 3; ++Index) DeltaP[Index] = T * Step[Index] + DeltaT * Ray[Index];
        Footprint[0] += std::fabs(Dot3(DeltaP, Ground.AxisU) * InvCell);
        Footprint[1] += std::fabs(Dot3(DeltaP, Ground.AxisV) * InvCell);
    }
    // 像素足迹（中心±半像素）的角点射线须全部朝向平面，否则足迹触及平行奇点，取两色均值
    const float HorizonMargin = NdotRay + 0.5f * (std::fabs(Dot3(N, Ground.RayStepX)) + std::fabs(Dot3(N, Ground.RayStepY)));
    float Coverage = 0.5f;
    if(HorizonMargin < 0.f && std::isfinite(U) && std::isfinite(V) && std::isfinite(Footprint[0]) && std::isfinite(Footprint[1]))
        Coverage = 0.5f * (1.f + FilteredSquareWave(U, Footprint[0]) * FilteredSquareWave(V, Footprint[1]));

    Out.Depth = T;
    Out.U = U;
    Out.V = V;
    Out.Coverage = Coverage;
    for(int Index = 0; Index < 3; ++Index)
        Out.Albedo[Index] = Ground.AlbedoOdd[Index] + (Ground.AlbedoEven[Index] - Ground.AlbedoOdd[Index]) * Coverage;
    return true;
}
