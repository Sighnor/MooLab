#version 460 core

// 光源立方体：对应CPU rasterizer中 case LIGHT（color = clampv(255 * light_radiance)）。
// 深度规则同pbr.frag：CPU只接受 depth > z_near 的像素，深度以观察空间距离写入（按2的幂缩放到[0, 1)）。

noperspective in float VViewDepth;

layout(location = 5) uniform vec3 ULightRadiance;
layout(location = 8) uniform float UZNear;
layout(location = 9) uniform float UDepthScale;

layout(location = 0) out vec4 OutColor;

// 按位判断NaN：isnan可能被编译器按“无NaN”假设优化掉
bool MooIsNan(float X)
{
    return (floatBitsToUint(X) & 0x7FFFFFFFu) > 0x7F800000u;
}

void main()
{
    // CPU：NaN深度的比较恒为假，不绘制
    if(MooIsNan(VViewDepth) || !(VViewDepth > UZNear))
    {
        discard;
    }
    gl_FragDepth = VViewDepth * UDepthScale;
    OutColor = vec4(clamp(ULightRadiance, 0.0, 1.0), 1.0);
}
