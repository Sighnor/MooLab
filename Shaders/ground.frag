#version 460 core

// 程序化棋盘地面：逐项对应CPU的EvaluateGroundPixel（Source/Runtime/Renderer/Private/PbrShading.cpp），
// 着色由拼接的 pbr_common.glsl 提供。与fullscreen.vert配合，在标准全幅视口、无scissor下绘制；
// 像素坐标由gl_FragCoord换算为图像列/行（行0为顶部），射线参数已在C++中包含模型光栅的一行上移，此处不再偏移。

in vec2 VUv;

layout(location = 0) uniform vec4 UPlane;
layout(location = 1) uniform vec3 UOrigin;
layout(location = 2) uniform vec3 UAxisU;
layout(location = 3) uniform vec3 UAxisV;
layout(location = 4) uniform vec3 URayTopLeft;
layout(location = 5) uniform vec3 URayStepX;
layout(location = 6) uniform vec3 URayStepY;
layout(location = 7) uniform float UCellSize;
layout(location = 8) uniform vec3 UAlbedoEven;
layout(location = 9) uniform vec3 UAlbedoOdd;
layout(location = 10) uniform float URoughness;
layout(location = 11) uniform float UMetallic;
layout(location = 12) uniform vec3 UCameraPos;
layout(location = 13) uniform vec3 ULightPos;
layout(location = 14) uniform vec3 ULightRadiance;
layout(location = 15) uniform float UZNear;
layout(location = 16) uniform float UClearDepth;
layout(location = 17) uniform float UDepthScale;
layout(location = 18) uniform float UTargetHeight;

layout(location = 0) out vec4 OutColor;

const float GroundMinFootprint = 1e-3;
const float GroundMaxFootprint = 64.0;

// 周期2的平方波（floor(x)偶为+1、奇为-1）的原函数
float MooSquareWaveIntegral(float X)
{
    return 1.0 - abs(X - 2.0 * floor(X / 2.0) - 1.0);
}

// 平方波在[X-W/2, X+W/2]上的盒滤波均值
float MooFilteredSquareWave(float X, float W)
{
    if(W >= GroundMaxFootprint)
    {
        return 0.0;
    }
    if(W < GroundMinFootprint)
    {
        W = GroundMinFootprint;
    }
    return MooClampf((MooSquareWaveIntegral(X + 0.5 * W) - MooSquareWaveIntegral(X - 0.5 * W)) / W, -1.0, 1.0);
}

bool MooIsFinite(float X)
{
    return (floatBitsToUint(X) & 0x7F800000u) != 0x7F800000u;
}

void main()
{
    vec3 N = UPlane.xyz;
    float Side = dot(N, UCameraPos) + UPlane.w;
    if(!(Side > 0.0))
    {
        discard;
    }

    float I = gl_FragCoord.x - 0.5;
    float J = UTargetHeight - 0.5 - gl_FragCoord.y;
    vec3 Ray = URayTopLeft + I * URayStepX + J * URayStepY;
    // 只接受朝向平面的射线；平行（含分母为0）与NaN一并拒绝
    float NdotRay = dot(N, Ray);
    if(!(NdotRay < 0.0))
    {
        discard;
    }
    float T = -Side / NdotRay;
    if(!(T > UZNear) || !(T <= UClearDepth))
    {
        discard;
    }

    vec3 Relative = (UCameraPos - UOrigin) + T * Ray;
    vec3 Position = UCameraPos + T * Ray;
    float InvCell = 1.0 / UCellSize;
    float U = dot(Relative, UAxisU) * InvCell;
    float V = dot(Relative, UAxisV) * InvCell;

    // 交点对像素轴的解析导数：dt = -t*(n·dq)/(n·q)，dP = t*dq + dt*q（不用fwidth，与CPU一致）
    vec2 Footprint = vec2(0.0);
    for(int Axis = 0; Axis < 2; Axis++)
    {
        vec3 Step = Axis == 0 ? URayStepX : URayStepY;
        float DeltaT = -T * dot(N, Step) / NdotRay;
        vec3 DeltaP = T * Step + DeltaT * Ray;
        Footprint.x += abs(dot(DeltaP, UAxisU) * InvCell);
        Footprint.y += abs(dot(DeltaP, UAxisV) * InvCell);
    }
    // 像素足迹（中心±半像素）的角点射线须全部朝向平面，否则足迹触及平行奇点，取两色均值
    float HorizonMargin = NdotRay + 0.5 * (abs(dot(N, URayStepX)) + abs(dot(N, URayStepY)));
    float Coverage = 0.5;
    if(HorizonMargin < 0.0 && MooIsFinite(U) && MooIsFinite(V) && MooIsFinite(Footprint.x) && MooIsFinite(Footprint.y))
    {
        Coverage = 0.5 * (1.0 + MooFilteredSquareWave(U, Footprint.x) * MooFilteredSquareWave(V, Footprint.y));
    }

    vec3 Albedo = UAlbedoOdd + (UAlbedoEven - UAlbedoOdd) * Coverage;
    vec3 Color = MooShadePbr(Albedo, N, Position, UCameraPos, ULightPos, ULightRadiance, URoughness, UMetallic);
    gl_FragDepth = T * UDepthScale;
    OutColor = vec4(MooNanToZero(Color), 1.0);
}
