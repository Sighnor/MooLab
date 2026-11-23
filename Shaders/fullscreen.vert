#version 460 core

// 覆盖全屏的单个三角形（无顶点缓冲，按gl_VertexID生成），VUv在屏幕内取值[0, 1]，原点在左下

out vec2 VUv;

void main()
{
    vec2 P = vec2(float((gl_VertexID << 1) & 2), float(gl_VertexID & 2));
    VUv = P;
    gl_Position = vec4(P * 2.0 - 1.0, 0.0, 1.0);
}
