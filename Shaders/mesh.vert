#version 460 core

// 网格顶点着色器：对应CPU rasterizer的三角形建立阶段（modelview = view * model，再乘投影）。
// 投影后的除法与CPU standardize相同（除以w + 1e-8），并输出w = 1：GL不再做透视除法与基于w的裁剪，
// w <= 0（相机后方）顶点也按CPU的方式“镜像”投影，覆盖区域与CPU一致；近平面由片元按深度丢弃。
// 属性以noperspective传给片元：CPU在屏幕空间按重心坐标线性插值（非透视校正），此处保持一致；
// 片元阶段再用model矩阵变换插值后的模型空间位置/法线，与CPU get_fragment_payload的运算次序相同。
// CPU深度为观察空间距离-z在屏幕空间线性插值的结果（非GL的NDC深度），一并以noperspective传出，由片元写入gl_FragDepth。

layout(location = 0) in vec3 InPosition;
layout(location = 1) in vec3 InNormal;
layout(location = 2) in vec2 InTexCoord;

layout(location = 0) uniform mat4 UModel;
layout(location = 1) uniform mat4 UView;
layout(location = 2) uniform mat4 UProj;

noperspective out vec3 VLocalPos;
noperspective out vec3 VLocalNormal;
noperspective out vec2 VTexCoord;
noperspective out float VViewDepth;

void main()
{
    mat4 ModelView = UView * UModel;
    vec4 ViewPos = ModelView * vec4(InPosition, 1.0);
    vec4 ClipPos = UProj * ViewPos;
    gl_Position = vec4(ClipPos.xy / (ClipPos.w + 1e-8), 0.0, 1.0);

    VLocalPos = InPosition;
    VLocalNormal = InNormal;
    VTexCoord = InTexCoord;
    VViewDepth = -ViewPos.z;
}
