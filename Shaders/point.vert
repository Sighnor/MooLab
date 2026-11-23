#version 460 core

// 像素方块：对应CPU draw_point。像素坐标由CPU按draw_point相同的公式算出（行0为图像顶部），
// 此处只把点放到该像素中心，以边长 2 * HalfSize + 1 的点精灵覆盖与CPU循环完全相同的像素集合。

layout(location = 0) in ivec2 InPixel;

layout(location = 0) uniform vec2 UTargetSize;
layout(location = 1) uniform float UPointSize;

void main()
{
    // 图像行0为顶部，GL窗口坐标行0为底部
    vec2 WindowPos = vec2(float(InPixel.x) + 0.5, UTargetSize.y - 1.0 - float(InPixel.y) + 0.5);
    gl_Position = vec4(WindowPos / UTargetSize * 2.0 - 1.0, 0.0, 1.0);
    gl_PointSize = UPointSize;
}
