#version 460 core

// PBR片元着色器：对应 shader.hpp 的 MooShader::pbr_shader。BRDF与输出变换由拼接的 pbr_common.glsl 提供
// （标准单次散射 metallic-roughness，与CPU数值核同一公式），本文件只保留网格插值、深度规则与存量纹理寻址
// （vec.hpp normalize带1e-8、texture::getcolor的行列寻址与float_to_int取整）。

noperspective in vec3 VLocalPos;
noperspective in vec3 VLocalNormal;
noperspective in vec2 VTexCoord;
noperspective in float VViewDepth;

layout(location = 0) uniform mat4 UModel;
layout(location = 3) uniform vec3 UCameraPos;
layout(location = 4) uniform vec3 ULightPos;
layout(location = 5) uniform vec3 ULightRadiance;
layout(location = 6) uniform float URoughness;
layout(location = 7) uniform float UMetallic;
layout(location = 8) uniform float UZNear;
// CPU深度缓冲存观察空间距离（清屏值100）；乘以2的幂缩放到[0, 1]写入D32F，缩放无舍入，比较结果与CPU一致
layout(location = 9) uniform float UDepthScale;

// 纹理以RGB32F存储CPU原值（0~255）
layout(binding = 0) uniform sampler2D UAlbedoMap;

layout(location = 0) out vec4 OutColor;

// global.hpp float_to_int：小数部分大于0.5才进位（0.5向下）
int MooFloatToInt(float F)
{
    float A = floor(F);
    float B = ceil(F);
    if(F - A > 0.5)
    {
        return int(B);
    }
    return int(A);
}

// texture::getcolor(i, j)：第一个参数乘行数取行，第二个参数乘列数取列（即UV互换），最近邻，返回0~255
vec3 MooSampleLegacy(sampler2D Tex, float I, float J)
{
    ivec2 Size = textureSize(Tex, 0);
    float Rows = float(Size.y);
    float Cols = float(Size.x);
    int Y = MooFloatToInt(MooClampf(I * Rows, 0.0, Rows - 1.0));
    int X = MooFloatToInt(MooClampf(J * Cols, 0.0, Cols - 1.0));
    return texelFetch(Tex, ivec2(X, Y), 0).rgb;
}

void main()
{
    // rasterizer：depth <= fbo深度（由GL_LEQUAL完成）且 depth > z_near；NaN深度的比较恒为假，不绘制
    if(MooIsNan(VViewDepth) || !(VViewDepth > UZNear))
    {
        discard;
    }
    gl_FragDepth = VViewDepth * UDepthScale;

    // get_fragment_payload：插值后的模型空间量经model矩阵变换，再按standardize（除以w + 1e-8）
    vec4 WorldPos4 = UModel * vec4(VLocalPos, 1.0);
    vec3 FragPos = (WorldPos4 / (WorldPos4.w + 1e-8)).xyz;
    vec4 WorldNormal4 = UModel * vec4(VLocalNormal, 0.0);
    vec3 FragNormal = (WorldNormal4 / (WorldNormal4.w + 1e-8)).xyz;

    // 纹理颜色按gamma 2.2解码为线性albedo
    vec3 Albedo = MooPowGamma(MooSampleLegacy(UAlbedoMap, VTexCoord.x, VTexCoord.y) / 255.0);
    vec3 N = MooNormalize(FragNormal);
    vec3 Color = MooShadePbr(Albedo, N, FragPos, UCameraPos, ULightPos, ULightRadiance, URoughness, UMetallic);

    // CPU：clampv(255 * color)后转8位；UNORM写出同样截断到[0, 1]
    OutColor = vec4(MooNanToZero(Color), 1.0);
}
