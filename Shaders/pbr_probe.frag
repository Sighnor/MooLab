#version 460 core

// material数值探针：逐列读取一组BRDF输入，调用拼接的 pbr_common.glsl 中与生产网格/地面相同的
// MooEvaluatePbrBrdf，把漫反射与镜面项原样写入两个RGBA32F目标（不做tonemap/量化）。仅供MaterialProbe使用。
// 输入纹理宽为样本数、高为3：行0 = (Albedo.rgb, Roughness)，行1 = (Metallic, NdotV, NdotL, NdotH)，行2 = (VdotH, -, -, -)。

layout(binding = 0) uniform sampler2D UInputs;

layout(location = 0) out vec4 OutDiffuse;
layout(location = 1) out vec4 OutSpecular;

void main()
{
    int Index = int(gl_FragCoord.x);
    vec4 A = texelFetch(UInputs, ivec2(Index, 0), 0);
    vec4 B = texelFetch(UInputs, ivec2(Index, 1), 0);
    vec4 C = texelFetch(UInputs, ivec2(Index, 2), 0);
    vec3 Diffuse;
    vec3 Specular;
    MooEvaluatePbrBrdf(A.rgb, A.a, B.x, B.y, B.z, B.w, C.x, Diffuse, Specular);
    OutDiffuse = vec4(Diffuse, 1.0);
    OutSpecular = vec4(Specular, 1.0);
}
