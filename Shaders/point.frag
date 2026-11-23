#version 460 core

layout(location = 2) uniform vec3 UColor;

layout(location = 0) out vec4 OutColor;

void main()
{
    OutColor = vec4(UColor, 1.0);
}
