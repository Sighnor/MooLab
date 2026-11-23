#version 460 core

// CPU参考路径呈现：纹理为CPU帧缓冲原值（RGB 0~255，行0为图像顶部），翻转行序后归一化输出

in vec2 VUv;

layout(binding = 0) uniform sampler2D UImage;

layout(location = 0) out vec4 OutColor;

void main()
{
    vec3 Color = texture(UImage, vec2(VUv.x, 1.0 - VUv.y)).rgb / 255.0;
    OutColor = vec4(clamp(Color, 0.0, 1.0), 1.0);
}
