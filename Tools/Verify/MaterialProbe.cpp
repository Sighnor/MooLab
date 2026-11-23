#include "AppWindow.h"
#include "GlRhi.h"
#include "PbrShading.h"
#include "Core.h"
#include <cmath>
#include <cstdio>
#include <limits>
#include <random>
#include <string>
#include <vector>

// material数值探针（验收工具，不进入生产目标）：同一组BRDF输入分别经CPU数值核EvaluatePbrBrdf与
// 拼接pbr_common.glsl的GPU探针着色器求值，GPU写RGBA32F并原样回读（不tonemap、不量化）。
// stdout逐样本输出 SAMPLE 行，由 verify_gl.py material 以Python double独立求值比对；
// 随机互易/有限/非负、背面与NaN、越界钳制等不依赖独立求值器的性质在此直接检查，失败即非零退出。

namespace
{
struct FProbeSample
{
    std::string Tag;
    FPbrBrdfInput Input;
};

FPbrBrdfInput MakeInput(float R, float G, float B, float Roughness, float Metallic,
                        float NdotV, float NdotL, float NdotH, float VdotH)
{
    FPbrBrdfInput Input;
    Input.Albedo[0] = R;
    Input.Albedo[1] = G;
    Input.Albedo[2] = B;
    Input.Roughness = Roughness;
    Input.Metallic = Metallic;
    Input.NdotV = NdotV;
    Input.NdotL = NdotL;
    Input.NdotH = NdotH;
    Input.VdotH = VdotH;
    return Input;
}

// 由法线(0,0,1)上方的V、L（double）求四个余弦；H = normalize(V+L)
void SetDirections(FPbrBrdfInput& Input, const double V[3], const double L[3])
{
    const double VLength = std::sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
    const double LLength = std::sqrt(L[0] * L[0] + L[1] * L[1] + L[2] * L[2]);
    double Vn[3], Ln[3], H[3];
    for(int Index = 0; Index < 3; ++Index)
    {
        Vn[Index] = V[Index] / VLength;
        Ln[Index] = L[Index] / LLength;
        H[Index] = Vn[Index] + Ln[Index];
    }
    const double HLength = std::sqrt(H[0] * H[0] + H[1] * H[1] + H[2] * H[2]);
    for(double& Value : H) Value /= HLength;
    Input.NdotV = static_cast<float>(Vn[2]);
    Input.NdotL = static_cast<float>(Ln[2]);
    Input.NdotH = static_cast<float>(H[2]);
    Input.VdotH = static_cast<float>(Vn[0] * H[0] + Vn[1] * H[1] + Vn[2] * H[2]);
}

void SphericalDirection(double ThetaDegrees, double PhiDegrees, double Out[3])
{
    const double Degree = 3.14159265358979323846 / 180.;
    Out[0] = std::sin(ThetaDegrees * Degree) * std::cos(PhiDegrees * Degree);
    Out[1] = std::sin(ThetaDegrees * Degree) * std::sin(PhiDegrees * Degree);
    Out[2] = std::cos(ThetaDegrees * Degree);
}

bool IsZero(const FPbrBrdfOutput& Output)
{
    for(int Channel = 0; Channel < 3; ++Channel)
        if(Output.Diffuse[Channel] != 0.f || Output.Specular[Channel] != 0.f) return false;
    return true;
}

bool SameOutput(const FPbrBrdfOutput& A, const FPbrBrdfOutput& B)
{
    for(int Channel = 0; Channel < 3; ++Channel)
        if(A.Diffuse[Channel] != B.Diffuse[Channel] || A.Specular[Channel] != B.Specular[Channel]) return false;
    return true;
}

std::vector<FProbeSample> MakeSamples()
{
    std::vector<FProbeSample> Samples;
    const float Cos30 = static_cast<float>(std::cos(3.14159265358979323846 / 6.));
    // 数值锚点：正入射三组 + L=N、NdotV=0.5（区分VdotH与NdotV的Fresnel）
    Samples.push_back({"anchor_dielectric", MakeInput(1.f, 1.f, 1.f, 0.5f, 0.f, 1.f, 1.f, 1.f, 1.f)});
    Samples.push_back({"anchor_metal", MakeInput(1.f, 1.f, 1.f, 0.5f, 1.f, 1.f, 1.f, 1.f, 1.f)});
    Samples.push_back({"anchor_half", MakeInput(0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 1.f, 1.f, 1.f, 1.f)});
    Samples.push_back({"anchor_fresnel", MakeInput(0.f, 0.f, 0.f, 1.f, 0.f, 0.5f, 1.f, Cos30, Cos30)});
    // 网格：m × r × 四组有效方向（正入射、45°镜面、离镜面、掠射）；彩色albedo
    const double Directions[4][4] = {{0., 0., 0., 0.}, {45., 0., 45., 180.}, {30., 0., 60., 90.}, {80., 0., 70., 150.}};
    for(float Metallic : {0.f, 0.5f, 1.f})
    {
        for(float Roughness : {0.05f, 0.5f, 1.f})
        {
            for(int Direction = 0; Direction < 4; ++Direction)
            {
                FProbeSample Sample;
                Sample.Tag = "grid_m" + std::to_string(static_cast<int>(Metallic * 10.f)) + "_r"
                    + std::to_string(static_cast<int>(Roughness * 100.f + 0.5f)) + "_d" + std::to_string(Direction);
                Sample.Input = MakeInput(0.9f, 0.5f, 0.1f, Roughness, Metallic, 0.f, 0.f, 0.f, 0.f);
                double V[3], L[3];
                SphericalDirection(Directions[Direction][0], Directions[Direction][1], V);
                SphericalDirection(Directions[Direction][2], Directions[Direction][3], L);
                SetDirections(Sample.Input, V, L);
                Samples.push_back(Sample);
            }
        }
    }
    // 边界材质与零输出样本：越界值由数值核钳制；背面/掠射为零
    Samples.push_back({"bound_r005_m0", MakeInput(0.3f, 0.6f, 0.9f, 0.05f, 0.f, 0.7f, 0.6f, 0.95f, 0.9f)});
    Samples.push_back({"bound_r1_m1", MakeInput(0.3f, 0.6f, 0.9f, 1.f, 1.f, 0.7f, 0.6f, 0.95f, 0.9f)});
    Samples.push_back({"black_dielectric", MakeInput(0.f, 0.f, 0.f, 0.3f, 0.f, 0.8f, 0.7f, 0.98f, 0.9f)});
    Samples.push_back({"zero_backface_view", MakeInput(0.5f, 0.5f, 0.5f, 0.5f, 0.f, -0.5f, 0.5f, 0.5f, 0.5f)});
    Samples.push_back({"zero_backface_light", MakeInput(0.5f, 0.5f, 0.5f, 0.5f, 0.f, 0.5f, -0.5f, 0.5f, 0.5f)});
    Samples.push_back({"zero_grazing", MakeInput(0.5f, 0.5f, 0.5f, 0.5f, 0.f, 0.f, 0.5f, 0.5f, 0.5f)});
    return Samples;
}

// GPU求值：输入纹理宽为样本数、高为3，两个RGBA32F目标各一行
void EvaluateGpu(const std::vector<FPbrBrdfInput>& Inputs, std::vector<FPbrBrdfOutput>& Outputs)
{
    const int Count = static_cast<int>(Inputs.size());
    std::vector<float> Texels(static_cast<size_t>(Count) * 3 * 4, 0.f);
    for(int Index = 0; Index < Count; ++Index)
    {
        const FPbrBrdfInput& Input = Inputs[Index];
        float* Row0 = &Texels[(static_cast<size_t>(0) * Count + Index) * 4];
        float* Row1 = &Texels[(static_cast<size_t>(1) * Count + Index) * 4];
        float* Row2 = &Texels[(static_cast<size_t>(2) * Count + Index) * 4];
        Row0[0] = Input.Albedo[0];
        Row0[1] = Input.Albedo[1];
        Row0[2] = Input.Albedo[2];
        Row0[3] = Input.Roughness;
        Row1[0] = Input.Metallic;
        Row1[1] = Input.NdotV;
        Row1[2] = Input.NdotL;
        Row1[3] = Input.NdotH;
        Row2[0] = Input.VdotH;
    }

    FGlProgram Program = CreateGlProgramFromFiles(GetGlShaderPath("fullscreen.vert"), GetGlShaderPath("pbr_probe.frag"),
                                                  GetGlShaderPath("pbr_common.glsl"));
    FGlVertexArray VertexArray = CreateGlVertexArray();
    FGlTexture2D InputTexture = CreateGlTexture2D(GL_RGBA32F, Count, 3);
    FGlTexture2D DiffuseTexture = CreateGlTexture2D(GL_RGBA32F, Count, 1);
    FGlTexture2D SpecularTexture = CreateGlTexture2D(GL_RGBA32F, Count, 1);
    FGlFramebuffer Framebuffer = CreateGlFramebuffer();

    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
    glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
    glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
    glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
    glTextureSubImage2D(InputTexture.Get(), 0, 0, 0, Count, 3, GL_RGBA, GL_FLOAT, Texels.data());

    glNamedFramebufferTexture(Framebuffer.Get(), GL_COLOR_ATTACHMENT0, DiffuseTexture.Get(), 0);
    glNamedFramebufferTexture(Framebuffer.Get(), GL_COLOR_ATTACHMENT1, SpecularTexture.Get(), 0);
    const GLenum DrawBuffers[] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1};
    glNamedFramebufferDrawBuffers(Framebuffer.Get(), 2, DrawBuffers);
    MOO_CHECK_ALWAYS(glCheckNamedFramebufferStatus(Framebuffer.Get(), GL_DRAW_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE,
                     "material probe framebuffer incomplete");

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, Framebuffer.Get());
    glViewport(0, 0, Count, 1);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_SCISSOR_TEST);
    glDisable(GL_BLEND);
    glDisable(GL_CULL_FACE);
    glUseProgram(Program.Get());
    glBindTextureUnit(0, InputTexture.Get());
    glBindVertexArray(VertexArray.Get());
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);

    std::vector<float> Diffuse(static_cast<size_t>(Count) * 4), Specular(static_cast<size_t>(Count) * 4);
    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
    glPixelStorei(GL_PACK_ALIGNMENT, 4);
    glPixelStorei(GL_PACK_ROW_LENGTH, 0);
    glPixelStorei(GL_PACK_SKIP_ROWS, 0);
    glPixelStorei(GL_PACK_SKIP_PIXELS, 0);
    glGetTextureImage(DiffuseTexture.Get(), 0, GL_RGBA, GL_FLOAT, static_cast<GLsizei>(Diffuse.size() * sizeof(float)), Diffuse.data());
    glGetTextureImage(SpecularTexture.Get(), 0, GL_RGBA, GL_FLOAT, static_cast<GLsizei>(Specular.size() * sizeof(float)), Specular.data());
    MOO_CHECK_ALWAYS(glGetError() == GL_NO_ERROR, "GL error in material probe");

    Outputs.assign(Inputs.size(), FPbrBrdfOutput());
    for(int Index = 0; Index < Count; ++Index)
    {
        for(int Channel = 0; Channel < 3; ++Channel)
        {
            Outputs[Index].Diffuse[Channel] = Diffuse[static_cast<size_t>(Index) * 4 + Channel];
            Outputs[Index].Specular[Channel] = Specular[static_cast<size_t>(Index) * 4 + Channel];
        }
    }
}

// 不依赖独立求值器的CPU性质：随机方向互易/有限/非负、边界钳制、零输出与金属漫反射
void CheckCpuProperties()
{
    std::mt19937 Random(20240607u);
    std::uniform_real_distribution<double> Unit(0., 1.);
    const double TwoPi = 2. * 3.14159265358979323846;
    for(int Iteration = 0; Iteration < 10000; ++Iteration)
    {
        // 上半球均匀方向（z ∈ (0,1]）
        double V[3], L[3];
        for(double* Direction : {V, L})
        {
            const double Z = 1. - Unit(Random) * 0.999;
            const double Phi = Unit(Random) * TwoPi;
            const double Radius = std::sqrt(1. - Z * Z);
            Direction[0] = Radius * std::cos(Phi);
            Direction[1] = Radius * std::sin(Phi);
            Direction[2] = Z;
        }
        FPbrBrdfInput Input = MakeInput(static_cast<float>(Unit(Random)), static_cast<float>(Unit(Random)), static_cast<float>(Unit(Random)),
                                        static_cast<float>(0.05 + 0.95 * Unit(Random)), static_cast<float>(Unit(Random)), 0.f, 0.f, 0.f, 0.f);
        SetDirections(Input, V, L);
        FPbrBrdfInput Swapped = Input;
        Swapped.NdotV = Input.NdotL;
        Swapped.NdotL = Input.NdotV;
        const FPbrBrdfOutput A = EvaluatePbrBrdf(Input);
        const FPbrBrdfOutput B = EvaluatePbrBrdf(Swapped);
        for(int Channel = 0; Channel < 3; ++Channel)
        {
            for(float Value : {A.Diffuse[Channel], A.Specular[Channel]})
                MOO_CHECK_ALWAYS(std::isfinite(Value) && Value >= 0.f, "random sample %d not finite/non-negative", Iteration);
            MOO_CHECK_ALWAYS(A.Diffuse[Channel] == B.Diffuse[Channel], "random sample %d diffuse not reciprocal", Iteration);
            const float Scale = std::fabs(A.Specular[Channel]) + std::fabs(B.Specular[Channel]);
            MOO_CHECK_ALWAYS(std::fabs(A.Specular[Channel] - B.Specular[Channel]) <= 1e-6f * Scale,
                             "random sample %d specular not reciprocal", Iteration);
        }
    }
    std::puts("cpu random 10000 reciprocity/finite/non-negative: PASS");

    const float Nan = std::numeric_limits<float>::quiet_NaN();
    const FPbrBrdfInput Base = MakeInput(0.6f, 0.4f, 0.2f, 0.4f, 0.3f, 0.7f, 0.6f, 0.95f, 0.9f);
    for(float FPbrBrdfInput::*Cosine : {&FPbrBrdfInput::NdotV, &FPbrBrdfInput::NdotL})
    {
        for(float Value : {0.f, -0.f, -0.5f, -1.f, Nan})
        {
            FPbrBrdfInput Input = Base;
            Input.*Cosine = Value;
            MOO_CHECK_ALWAYS(IsZero(EvaluatePbrBrdf(Input)), "back-facing/zero/NaN cosine did not return zero");
        }
    }
    std::puts("cpu back-face / zero-vector / NaN cosine zero: PASS");

    FPbrBrdfInput Low = Base, LowClamped = Base, High = Base, HighClamped = Base;
    Low.Roughness = 0.f;
    LowClamped.Roughness = 0.05f;
    Low.Metallic = -1.f;
    LowClamped.Metallic = 0.f;
    High.Roughness = 3.f;
    HighClamped.Roughness = 1.f;
    High.Metallic = 2.f;
    HighClamped.Metallic = 1.f;
    MOO_CHECK_ALWAYS(SameOutput(EvaluatePbrBrdf(Low), EvaluatePbrBrdf(LowClamped)), "lower material clamp wrong");
    MOO_CHECK_ALWAYS(SameOutput(EvaluatePbrBrdf(High), EvaluatePbrBrdf(HighClamped)), "upper material clamp wrong");
    FPbrBrdfInput OverCosine = Base, UnitCosine = Base;
    OverCosine.NdotH = 1.5f;
    UnitCosine.NdotH = 1.f;
    MOO_CHECK_ALWAYS(SameOutput(EvaluatePbrBrdf(OverCosine), EvaluatePbrBrdf(UnitCosine)), "cosine clamp wrong");
    std::puts("cpu roughness/metallic/cosine clamp: PASS");

    FPbrBrdfInput Metal = Base;
    Metal.Metallic = 1.f;
    const FPbrBrdfOutput MetalOutput = EvaluatePbrBrdf(Metal);
    FPbrBrdfInput Black = Base;
    Black.Metallic = 0.f;
    Black.Albedo[0] = Black.Albedo[1] = Black.Albedo[2] = 0.f;
    const FPbrBrdfOutput BlackOutput = EvaluatePbrBrdf(Black);
    for(int Channel = 0; Channel < 3; ++Channel)
    {
        MOO_CHECK_ALWAYS(MetalOutput.Diffuse[Channel] == 0.f, "metallic=1 diffuse not exactly zero");
        MOO_CHECK_ALWAYS(BlackOutput.Diffuse[Channel] == 0.f && BlackOutput.Specular[Channel] > 0.f, "black dielectric specular vanished");
    }
    std::puts("cpu metal diffuse zero / black dielectric specular: PASS");
}
}

int main()
{
    CheckCpuProperties();

    FAppWindow Window({64, 64, "MOOLAB material probe", true});
    CheckGlContextReady(Window);
    const std::vector<FProbeSample> Samples = MakeSamples();
    std::vector<FPbrBrdfInput> Inputs;
    for(const FProbeSample& Sample : Samples) Inputs.push_back(Sample.Input);
    // GPU的NaN余弦同样须走零分支（只在此检查，不输出到SAMPLE行）
    const float Nan = std::numeric_limits<float>::quiet_NaN();
    Inputs.push_back(MakeInput(0.5f, 0.5f, 0.5f, 0.5f, 0.f, Nan, 0.5f, 0.5f, 0.5f));
    Inputs.push_back(MakeInput(0.5f, 0.5f, 0.5f, 0.5f, 0.f, 0.5f, Nan, 0.5f, 0.5f));
    std::vector<FPbrBrdfOutput> GpuOutputs;
    EvaluateGpu(Inputs, GpuOutputs);
    for(size_t Index = Samples.size(); Index < Inputs.size(); ++Index)
        MOO_CHECK_ALWAYS(IsZero(GpuOutputs[Index]), "GPU NaN cosine did not return zero");
    std::puts("gpu NaN cosine zero: PASS");

    for(size_t Index = 0; Index < Samples.size(); ++Index)
    {
        const FPbrBrdfInput& In = Samples[Index].Input;
        const FPbrBrdfOutput Cpu = EvaluatePbrBrdf(In);
        const FPbrBrdfOutput& Gpu = GpuOutputs[Index];
        for(int Channel = 0; Channel < 3; ++Channel)
        {
            for(float Value : {Gpu.Diffuse[Channel], Gpu.Specular[Channel]})
                MOO_CHECK_ALWAYS(std::isfinite(Value), "GPU output not finite for %s", Samples[Index].Tag.c_str());
        }
        std::printf("SAMPLE %s %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g"
                    " %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g\n",
                    Samples[Index].Tag.c_str(), In.Albedo[0], In.Albedo[1], In.Albedo[2], In.Roughness, In.Metallic,
                    In.NdotV, In.NdotL, In.NdotH, In.VdotH,
                    Cpu.Diffuse[0], Cpu.Diffuse[1], Cpu.Diffuse[2], Cpu.Specular[0], Cpu.Specular[1], Cpu.Specular[2],
                    Gpu.Diffuse[0], Gpu.Diffuse[1], Gpu.Diffuse[2], Gpu.Specular[0], Gpu.Specular[1], Gpu.Specular[2]);
    }
    std::printf("material-probe samples=%zu: PASS\n", Samples.size());
    return 0;
}
