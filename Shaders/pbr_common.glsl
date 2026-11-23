// 片元公共源码：由CreateGlProgramFromFiles拼接在主文件版本指令之后，本文件不得含版本指令或uniform。
// 标准单次散射 metallic-roughness 数值核，与CPU的EvaluatePbrBrdf（Source/Runtime/Renderer/Private/PbrShading.cpp）
// 及MooShadePbr（Source/Legacy/shader.hpp）逐项同一公式与运算次序；修改任一处须同步另外两处并重跑material探针。

const float PI = 3.14159265358979323846;

// global.hpp clampf：用比较实现（与std::clamp/GLSL clamp对NaN的处理不同，保持CPU语义）
float MooClampf(float X, float MinValue, float MaxValue)
{
    if(X < MinValue)
    {
        X = MinValue;
    }
    if(X > MaxValue)
    {
        X = MaxValue;
    }
    return X;
}

// vec.hpp normalize：v / (length(v) + 1e-8)
vec3 MooNormalize(vec3 V)
{
    return V / (length(V) + 1e-8);
}

// 五次幂用连乘：GLSL的pow对负底数未定义
float MooPow5(float X)
{
    float X2 = X * X;
    return X2 * X2 * X;
}

// std::pow(x, 2.2f)：GLSL的pow按exp2(y * log2(x))定义精度，x = 0时log2未定义；CPU结果为0
vec3 MooPowGamma(vec3 X)
{
    return vec3(X.r > 0.0 ? pow(X.r, 2.2) : 0.0,
                X.g > 0.0 ? pow(X.g, 2.2) : 0.0,
                X.b > 0.0 ? pow(X.b, 2.2) : 0.0);
}

// 按位判断NaN：isnan可能被编译器按“无NaN”假设优化掉
bool MooIsNan(float X)
{
    return (floatBitsToUint(X) & 0x7FFFFFFFu) > 0x7F800000u;
}

// CPU的NaN经clampv保持NaN，再由OpenCV convertTo（cvRound后饱和）得到0；GL对NaN转UNORM未作规定，显式置0
vec3 MooNanToZero(vec3 C)
{
    return vec3(MooIsNan(C.r) ? 0.0 : C.r, MooIsNan(C.g) ? 0.0 : C.g, MooIsNan(C.b) ? 0.0 : C.b);
}

float MooSchlickG1(float X, float K)
{
    return X / (X * (1.0 - K) + K);
}

// 直接光BRDF：输出独立的漫反射与镜面项；Albedo为线性值。NdotV或NdotL非正（含NaN）时两项为零。
void MooEvaluatePbrBrdf(vec3 Albedo, float RoughnessIn, float MetallicIn,
                        float NdotVIn, float NdotLIn, float NdotHIn, float VdotHIn,
                        out vec3 Diffuse, out vec3 Specular)
{
    Diffuse = vec3(0.0);
    Specular = vec3(0.0);
    // NaN按位判断：驱动可能按“无NaN”假设把取反比较化简为 <= 0
    if(MooIsNan(NdotVIn) || MooIsNan(NdotLIn) || NdotVIn <= 0.0 || NdotLIn <= 0.0)
    {
        return;
    }

    float NdotV = MooClampf(NdotVIn, 0.0, 1.0);
    float NdotL = MooClampf(NdotLIn, 0.0, 1.0);
    float NdotH = MooClampf(NdotHIn, 0.0, 1.0);
    float VdotH = MooClampf(VdotHIn, 0.0, 1.0);
    float Roughness = MooClampf(RoughnessIn, 0.05, 1.0);
    float Metallic = MooClampf(MetallicIn, 0.0, 1.0);

    float Alpha = Roughness * Roughness;
    float Alpha2 = Alpha * Alpha;
    // precise禁止编译器把稳定形式重结合为 1 - NdotH2*(1-Alpha2)（低粗糙度、NdotH≈1时严重相消）
    precise float NdotH2 = NdotH * NdotH;
    precise float B = (1.0 - NdotH2) + NdotH2 * Alpha2;
    float D = Alpha2 / (PI * B * B);

    float RoughnessPlusOne = Roughness + 1.0;
    float K = RoughnessPlusOne * RoughnessPlusOne / 8.0;
    float G = MooSchlickG1(NdotV, K) * MooSchlickG1(NdotL, K);

    float X5 = MooPow5(1.0 - VdotH);
    float Denominator = 4.0 * NdotL * NdotV > 0.001 ? 4.0 * NdotL * NdotV : 0.001;

    vec3 F0 = (1.0 - Metallic) * 0.04 + Metallic * Albedo;
    vec3 F = F0 + (vec3(1.0) - F0) * X5;
    Specular = D * G * F / Denominator;
    Diffuse = (1.0 - Metallic) * (vec3(1.0) - F) * Albedo / PI;
}

// 点光着色（无距离衰减）：常量底光0.1 + BRDF*radiance*NdotL，再Reinhard与gamma 1/2.2；N须已归一化。
vec3 MooShadePbr(vec3 Albedo, vec3 N, vec3 FragPos, vec3 CameraPos, vec3 LightPos, vec3 LightRadiance,
                 float Roughness, float Metallic)
{
    vec3 V = MooNormalize(CameraPos - FragPos);
    vec3 L = MooNormalize(LightPos - FragPos);
    vec3 H = MooNormalize(V + L);

    vec3 Diffuse;
    vec3 Specular;
    MooEvaluatePbrBrdf(Albedo, Roughness, Metallic, dot(N, V), dot(N, L), dot(N, H), dot(V, H), Diffuse, Specular);
    vec3 BRDF = Diffuse + Specular;

    float NdotL = max(dot(N, L), 0.0);
    vec3 Color = vec3(0.1) + BRDF * LightRadiance * NdotL;

    Color = Color / (Color + vec3(1.0));
    Color = pow(Color, vec3(1.0 / 2.2));
    return Color;
}
