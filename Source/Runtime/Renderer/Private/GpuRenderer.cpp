#include "GpuRenderer.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

#include "Core.h"
#include "GlRhi.h"

namespace
{
// 着色器中以layout(location)固定的uniform位置，见shaders/*.vert|frag
constexpr GLint MeshUniform_Model = 0;
constexpr GLint MeshUniform_View = 1;
constexpr GLint MeshUniform_Proj = 2;
constexpr GLint MeshUniform_CameraPos = 3;
constexpr GLint MeshUniform_LightPos = 4;
constexpr GLint MeshUniform_LightRadiance = 5;
constexpr GLint MeshUniform_Roughness = 6;
constexpr GLint MeshUniform_Metallic = 7;
constexpr GLint MeshUniform_ZNear = 8;
constexpr GLint MeshUniform_DepthScale = 9;

constexpr GLint PointUniform_TargetSize = 0;
constexpr GLint PointUniform_PointSize = 1;
constexpr GLint PointUniform_Color = 2;

constexpr GLint GroundUniform_Plane = 0;
constexpr GLint GroundUniform_Origin = 1;
constexpr GLint GroundUniform_AxisU = 2;
constexpr GLint GroundUniform_AxisV = 3;
constexpr GLint GroundUniform_RayTopLeft = 4;
constexpr GLint GroundUniform_RayStepX = 5;
constexpr GLint GroundUniform_RayStepY = 6;
constexpr GLint GroundUniform_CellSize = 7;
constexpr GLint GroundUniform_AlbedoEven = 8;
constexpr GLint GroundUniform_AlbedoOdd = 9;
constexpr GLint GroundUniform_Roughness = 10;
constexpr GLint GroundUniform_Metallic = 11;
constexpr GLint GroundUniform_CameraPos = 12;
constexpr GLint GroundUniform_LightPos = 13;
constexpr GLint GroundUniform_LightRadiance = 14;
constexpr GLint GroundUniform_ZNear = 15;
constexpr GLint GroundUniform_ClearDepth = 16;
constexpr GLint GroundUniform_DepthScale = 17;
constexpr GLint GroundUniform_TargetHeight = 18;

constexpr GLuint TextureUnit_AlbedoMap = 0;

struct FGpuMesh
{
    FGlBuffer Positions;
    FGlBuffer Normals;
    FGlBuffer TexCoords;
    FGlBuffer Indices;
    FGlVertexArray VertexArray;
    int VertexCount = 0;
    GLsizei IndexCount = 0;
};

// 深度缩放：取2^-E使ClearDepth * Scale = m ∈ [0.5, 1)（ClearDepth = m * 2^E）。乘2的幂无舍入，比较结果与CPU一致；
// 清屏值严格小于1，被DEPTH_CLAMP夹到1.0的更远片元不会通过LEQUAL（ClearDepth恰为2的幂时也成立）
float GetDepthScale(float ClearDepth)
{
    int Exponent = 0;
    std::frexp(ClearDepth, &Exponent);
    return std::ldexp(1.f, -Exponent);
}
}

struct FGpuRenderer::FImpl
{
    int Width = 0;
    int Height = 0;

    FGlProgram LightProgram;
    FGlProgram PbrProgram;
    FGlProgram PointProgram;
    FGlProgram GroundProgram;

    std::vector<FGpuMesh> Meshes;
    std::vector<FGlTexture2D> Textures;

    FGlTexture2D ColorTarget;
    FGlTexture2D DepthTarget;
    FGlFramebuffer Framebuffer;

    FGlVertexArray PointVertexArray;
    FGlBuffer PointBuffer;
    int PointCapacity = 0;

    // 地面全屏三角形由gl_VertexID生成，空VAO仅满足core profile绘制要求
    FGlVertexArray GroundVertexArray;

    FGlGpuTimer Timer;

    FFrameParams Frame;
    float DepthScale = 1.f;
    float MaxPointSize = 1.f;
    bool bInFrame = false;
    // 地面须是帧内第一个绘制且至多一次：BeginFrame置位，任何绘制后清除
    bool bGroundAllowed = false;
};

FGpuRenderer::FGpuRenderer(const IGlContextOwner& Context, int Width, int Height)
{
    // FImpl中的计时器也会创建GL对象，必须在分配它之前验证上下文。
    CheckGlContextReady(Context);
    Impl = std::make_unique<FImpl>();
    MOO_CHECK_ALWAYS(Width > 0 && Height > 0, "invalid render target size %dx%d", Width, Height);
    Impl->Width = Width;
    Impl->Height = Height;

    Impl->LightProgram = CreateGlProgramFromFiles(GetGlShaderPath("mesh.vert"), GetGlShaderPath("light.frag"));
    Impl->PbrProgram = CreateGlProgramFromFiles(GetGlShaderPath("mesh.vert"), GetGlShaderPath("pbr.frag"), GetGlShaderPath("pbr_common.glsl"));
    Impl->PointProgram = CreateGlProgramFromFiles(GetGlShaderPath("point.vert"), GetGlShaderPath("point.frag"));
    Impl->GroundProgram = CreateGlProgramFromFiles(GetGlShaderPath("fullscreen.vert"), GetGlShaderPath("ground.frag"), GetGlShaderPath("pbr_common.glsl"));
    Impl->GroundVertexArray = CreateGlVertexArray();

    //离屏目标：截图与对拍只读这里，不受ImGui与窗口缩放影响
    Impl->ColorTarget = CreateGlTexture2D(GL_RGBA8, Width, Height);
    Impl->DepthTarget = CreateGlTexture2D(GL_DEPTH_COMPONENT32F, Width, Height);
    Impl->Framebuffer = CreateGlFramebuffer();
    glNamedFramebufferTexture(Impl->Framebuffer.Get(), GL_COLOR_ATTACHMENT0, Impl->ColorTarget.Get(), 0);
    glNamedFramebufferTexture(Impl->Framebuffer.Get(), GL_DEPTH_ATTACHMENT, Impl->DepthTarget.Get(), 0);
    const GLenum Status = glCheckNamedFramebufferStatus(Impl->Framebuffer.Get(), GL_FRAMEBUFFER);
    MOO_CHECK_ALWAYS(Status == GL_FRAMEBUFFER_COMPLETE, "offscreen framebuffer incomplete: 0x%x", Status);

    Impl->PointVertexArray = CreateGlVertexArray();
    glEnableVertexArrayAttrib(Impl->PointVertexArray.Get(), 0);
    glVertexArrayAttribIFormat(Impl->PointVertexArray.Get(), 0, 2, GL_INT, 0);
    glVertexArrayAttribBinding(Impl->PointVertexArray.Get(), 0, 0);

    GLfloat PointSizeRange[2] = {1.f, 1.f};
    glGetFloatv(GL_POINT_SIZE_RANGE, PointSizeRange);
    Impl->MaxPointSize = PointSizeRange[1];
}

FGpuRenderer::~FGpuRenderer() = default;

FMeshId FGpuRenderer::CreateMesh(const FMeshDataView& Mesh)
{
    MOO_CHECK_ALWAYS(Mesh.Positions && Mesh.Normals && Mesh.TexCoords && Mesh.Indices, "mesh view has null attribute pointers");
    MOO_CHECK_ALWAYS(Mesh.VertexCount > 0 && Mesh.TriangleCount > 0, "empty mesh (%d vertices, %d triangles)", Mesh.VertexCount, Mesh.TriangleCount);
    //越界索引会让GL读取缓冲外顶点（未定义行为），上传前一次性校验
    const int IndexCount = Mesh.TriangleCount * 3;
    const uint16_t* const IndexEnd = Mesh.Indices + IndexCount;
    const uint16_t* const BadIndex = std::find_if(Mesh.Indices, IndexEnd, [&](uint16_t Index) { return Index >= Mesh.VertexCount; });
    MOO_CHECK_ALWAYS(BadIndex == IndexEnd, "mesh index %d out of range (%d vertices)", BadIndex == IndexEnd ? 0 : *BadIndex, Mesh.VertexCount);

    const GLsizeiptr VertexCount = Mesh.VertexCount;
    FGpuMesh GpuMesh;
    GpuMesh.VertexCount = Mesh.VertexCount;
    GpuMesh.Positions = CreateGlBuffer(VertexCount * 3 * sizeof(float), Mesh.Positions, GL_DYNAMIC_STORAGE_BIT);
    GpuMesh.Normals = CreateGlBuffer(VertexCount * 3 * sizeof(float), Mesh.Normals, GL_DYNAMIC_STORAGE_BIT);
    GpuMesh.TexCoords = CreateGlBuffer(VertexCount * 2 * sizeof(float), Mesh.TexCoords, 0);
    GpuMesh.IndexCount = IndexCount;
    GpuMesh.Indices = CreateGlBuffer(static_cast<GLsizeiptr>(GpuMesh.IndexCount) * sizeof(uint16_t), Mesh.Indices, 0);

    //各属性独立缓冲（绑定点0/1/2对应location 0/1/2），后续compute阶段可单独改写位置/法线
    GpuMesh.VertexArray = CreateGlVertexArray();
    const GLuint Vao = GpuMesh.VertexArray.Get();
    glVertexArrayVertexBuffer(Vao, 0, GpuMesh.Positions.Get(), 0, 3 * sizeof(float));
    glVertexArrayVertexBuffer(Vao, 1, GpuMesh.Normals.Get(), 0, 3 * sizeof(float));
    glVertexArrayVertexBuffer(Vao, 2, GpuMesh.TexCoords.Get(), 0, 2 * sizeof(float));
    const GLint ComponentCounts[3] = {3, 3, 2};
    for(GLuint Attrib = 0; Attrib < 3; Attrib++)
    {
        glEnableVertexArrayAttrib(Vao, Attrib);
        glVertexArrayAttribFormat(Vao, Attrib, ComponentCounts[Attrib], GL_FLOAT, GL_FALSE, 0);
        glVertexArrayAttribBinding(Vao, Attrib, Attrib);
    }
    glVertexArrayElementBuffer(Vao, GpuMesh.Indices.Get());

    Impl->Meshes.push_back(std::move(GpuMesh));
    return FMeshId{static_cast<int>(Impl->Meshes.size()) - 1};
}

void FGpuRenderer::UpdateMeshVertices(FMeshId Mesh, const float* Positions, const float* Normals, int VertexCount)
{
    MOO_CHECK_ALWAYS(!Impl->bInFrame, "mesh updates must precede BeginFrame");
    MOO_CHECK_ALWAYS(Mesh.Index >= 0 && static_cast<size_t>(Mesh.Index) < Impl->Meshes.size(), "invalid mesh id %d", Mesh.Index);
    FGpuMesh& Target = Impl->Meshes[Mesh.Index];
    MOO_CHECK_ALWAYS(Positions && Normals && VertexCount == Target.VertexCount, "mesh update requires full position/normal arrays with unchanged vertex count");
    const GLsizeiptr Bytes = static_cast<GLsizeiptr>(VertexCount) * 3 * sizeof(float);
    glNamedBufferSubData(Target.Positions.Get(), 0, Bytes, Positions);
    glNamedBufferSubData(Target.Normals.Get(), 0, Bytes, Normals);
}

FTextureId FGpuRenderer::CreateTextureRGBF(const FImageRGBFView& Image)
{
    MOO_CHECK_ALWAYS(Image.Data && Image.Width > 0 && Image.Height > 0, "invalid texture image %dx%d", Image.Width, Image.Height);

    FGlTexture2D Texture = CreateGlTexture2D(GL_RGB32F, Image.Width, Image.Height);
    //CPU存储行0为图像首行，原样作为纹理行0上传；着色器按CPU的(行, 列)寻址规则texelFetch
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
    glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
    glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
    glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
    glTextureSubImage2D(Texture.Get(), 0, 0, 0, Image.Width, Image.Height, GL_RGB, GL_FLOAT, Image.Data);

    Impl->Textures.push_back(std::move(Texture));
    return FTextureId{static_cast<int>(Impl->Textures.size()) - 1};
}

void FGpuRenderer::BeginFrame(const FFrameParams& Params)
{
    MOO_CHECK(!Impl->bInFrame);
    MOO_CHECK_ALWAYS(std::isfinite(Params.ClearDepth) && Params.ClearDepth > 0.f, "clear depth must be finite and positive, got %f", Params.ClearDepth);
    Impl->bInFrame = true;
    Impl->bGroundAllowed = true;
    Impl->Frame = Params;
    Impl->DepthScale = GetDepthScale(Params.ClearDepth);

    Impl->Timer.Begin();

    glBindFramebuffer(GL_FRAMEBUFFER, Impl->Framebuffer.Get());
    glViewport(0, 0, Impl->Width, Impl->Height);
    //ImGui后端会改写并恢复部分状态，这里每帧显式设置本渲染器依赖的全部状态
    glDisable(GL_SCISSOR_TEST);
    glDisable(GL_BLEND);
    glDisable(GL_CULL_FACE);
    glDisable(GL_STENCIL_TEST);
    glDisable(GL_PRIMITIVE_RESTART);
    glDisable(GL_PRIMITIVE_RESTART_FIXED_INDEX);
    glDisable(GL_RASTERIZER_DISCARD);
    //GL允许UNORM转换在相邻两个整数间抖动，关闭以保证确定的最近舍入
    glDisable(GL_DITHER);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glDepthRange(0.0, 1.0);
    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
    glDepthMask(GL_TRUE);

    const GLfloat ClearColor[4] = {Params.ClearColor[0], Params.ClearColor[1], Params.ClearColor[2], 1.f};
    const GLfloat ClearDepth = Params.ClearDepth * Impl->DepthScale;
    glClearNamedFramebufferfv(Impl->Framebuffer.Get(), GL_COLOR, 0, ClearColor);
    glClearNamedFramebufferfv(Impl->Framebuffer.Get(), GL_DEPTH, 0, &ClearDepth);
}

void FGpuRenderer::DrawMesh(FMeshId Mesh, const float Model[16], EShading Shading, const FPbrParams* PbrParams)
{
    MOO_CHECK(Impl->bInFrame);
    MOO_CHECK_ALWAYS(Mesh.Index >= 0 && Mesh.Index < static_cast<int>(Impl->Meshes.size()), "invalid mesh id %d", Mesh.Index);
    Impl->bGroundAllowed = false;
    const FGpuMesh& GpuMesh = Impl->Meshes[Mesh.Index];
    const FFrameParams& Frame = Impl->Frame;

    //CPU深度规则：ZNear < depth <= 已存深度，后绘制者在相等时胜出；深度由片元写入，关闭近远裁剪以免剔除CPU会保留的片元
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_DEPTH_CLAMP);
    //复刻CPU光栅的像素映射：CPU采样行y（自下而上，NDC取行中心）写入图像行rows-2-y，比标准映射整体上移1行；
    //视口上移1行即得到相同映射。CPU包围盒上限clamp到cols-1/rows-1且循环用严格小于，
    //最后一列与（自下而上的）最后一采样行从不着色，故以scissor排除最后一列与图像底行
    glViewport(0, 1, Impl->Width, Impl->Height);
    glEnable(GL_SCISSOR_TEST);
    glScissor(0, 1, Impl->Width - 1, Impl->Height - 1);

    GLuint Program = 0;
    if(Shading == EShading::Pbr)
    {
        MOO_CHECK_ALWAYS(PbrParams != nullptr, "PBR shading requires material parameters");
        const int TextureCount = static_cast<int>(Impl->Textures.size());
        const FTextureId Texture = PbrParams->AlbedoMap;
        MOO_CHECK_ALWAYS(Texture.Index >= 0 && Texture.Index < TextureCount, "invalid texture id %d", Texture.Index);

        Program = Impl->PbrProgram.Get();
        glProgramUniform3fv(Program, MeshUniform_CameraPos, 1, Frame.CameraPos);
        glProgramUniform3fv(Program, MeshUniform_LightPos, 1, Frame.LightPos);
        glProgramUniform1f(Program, MeshUniform_Roughness, PbrParams->Roughness);
        glProgramUniform1f(Program, MeshUniform_Metallic, PbrParams->Metallic);
        glBindTextureUnit(TextureUnit_AlbedoMap, Impl->Textures[Texture.Index].Get());
    }
    else
    {
        Program = Impl->LightProgram.Get();
    }

    glProgramUniformMatrix4fv(Program, MeshUniform_Model, 1, GL_FALSE, Model);
    glProgramUniformMatrix4fv(Program, MeshUniform_View, 1, GL_FALSE, Frame.View);
    glProgramUniformMatrix4fv(Program, MeshUniform_Proj, 1, GL_FALSE, Frame.Proj);
    glProgramUniform3fv(Program, MeshUniform_LightRadiance, 1, Frame.LightRadiance);
    glProgramUniform1f(Program, MeshUniform_ZNear, Frame.ZNear);
    glProgramUniform1f(Program, MeshUniform_DepthScale, Impl->DepthScale);

    glUseProgram(Program);
    glBindVertexArray(GpuMesh.VertexArray.Get());
    glDrawElements(GL_TRIANGLES, GpuMesh.IndexCount, GL_UNSIGNED_SHORT, nullptr);
}

void FGpuRenderer::DrawPixelSquares(const int* PixelXY, int Count, int HalfSize, const float ColorRGB[3])
{
    MOO_CHECK(Impl->bInFrame);
    Impl->bGroundAllowed = false;
    if(Count <= 0)
    {
        return;
    }
    MOO_CHECK_ALWAYS(PixelXY != nullptr && HalfSize >= 0, "invalid pixel squares input");
    const float PointSize = static_cast<float>(2 * HalfSize + 1);
    MOO_CHECK_ALWAYS(PointSize <= Impl->MaxPointSize, "point size %.0f exceeds GL_POINT_SIZE_RANGE max %.0f", PointSize, Impl->MaxPointSize);

    if(Count > Impl->PointCapacity)
    {
        Impl->PointCapacity = Count * 2;
        Impl->PointBuffer = CreateGlBuffer(static_cast<GLsizeiptr>(Impl->PointCapacity) * 2 * sizeof(int), nullptr, GL_DYNAMIC_STORAGE_BIT);
        glVertexArrayVertexBuffer(Impl->PointVertexArray.Get(), 0, Impl->PointBuffer.Get(), 0, 2 * sizeof(int));
    }
    glNamedBufferSubData(Impl->PointBuffer.Get(), 0, static_cast<GLsizeiptr>(Count) * 2 * sizeof(int), PixelXY);

    //CPU draw_point直接覆盖颜色，不读写深度；像素坐标已是图像坐标，使用标准视口、不裁剪（可写最后一列与底行，同CPU）
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_PROGRAM_POINT_SIZE);
    glViewport(0, 0, Impl->Width, Impl->Height);
    glDisable(GL_SCISSOR_TEST);

    const GLuint Program = Impl->PointProgram.Get();
    glProgramUniform2f(Program, PointUniform_TargetSize, static_cast<float>(Impl->Width), static_cast<float>(Impl->Height));
    glProgramUniform1f(Program, PointUniform_PointSize, PointSize);
    glProgramUniform3fv(Program, PointUniform_Color, 1, ColorRGB);

    glUseProgram(Program);
    glBindVertexArray(Impl->PointVertexArray.Get());
    glDrawArrays(GL_POINTS, 0, Count);
}

void FGpuRenderer::DrawGround(const FGroundParams& Ground)
{
    MOO_CHECK_ALWAYS(Impl->bInFrame, "DrawGround outside BeginFrame/EndFrame");
    MOO_CHECK_ALWAYS(Impl->bGroundAllowed, "DrawGround must be the first draw of a frame and called at most once");
    MOO_CHECK_ALWAYS(Ground.TargetWidth == Impl->Width && Ground.TargetHeight == Impl->Height,
                     "ground dimensions %dx%d do not match target %dx%d", Ground.TargetWidth, Ground.TargetHeight, Impl->Width, Impl->Height);
    Impl->bGroundAllowed = false;
    const FFrameParams& Frame = Impl->Frame;

    //标准全幅视口、无scissor：片元由gl_FragCoord换算图像列/行，射线参数已含一行上移；可覆盖模型不绘制的最后一列与底行
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_DEPTH_CLAMP);
    glViewport(0, 0, Impl->Width, Impl->Height);
    glDisable(GL_SCISSOR_TEST);

    const GLuint Program = Impl->GroundProgram.Get();
    glProgramUniform4fv(Program, GroundUniform_Plane, 1, Ground.Plane);
    glProgramUniform3fv(Program, GroundUniform_Origin, 1, Ground.Origin);
    glProgramUniform3fv(Program, GroundUniform_AxisU, 1, Ground.AxisU);
    glProgramUniform3fv(Program, GroundUniform_AxisV, 1, Ground.AxisV);
    glProgramUniform3fv(Program, GroundUniform_RayTopLeft, 1, Ground.RayTopLeft);
    glProgramUniform3fv(Program, GroundUniform_RayStepX, 1, Ground.RayStepX);
    glProgramUniform3fv(Program, GroundUniform_RayStepY, 1, Ground.RayStepY);
    glProgramUniform1f(Program, GroundUniform_CellSize, Ground.CellSize);
    glProgramUniform3fv(Program, GroundUniform_AlbedoEven, 1, Ground.AlbedoEven);
    glProgramUniform3fv(Program, GroundUniform_AlbedoOdd, 1, Ground.AlbedoOdd);
    glProgramUniform1f(Program, GroundUniform_Roughness, Ground.Roughness);
    glProgramUniform1f(Program, GroundUniform_Metallic, Ground.Metallic);
    glProgramUniform3fv(Program, GroundUniform_CameraPos, 1, Frame.CameraPos);
    glProgramUniform3fv(Program, GroundUniform_LightPos, 1, Frame.LightPos);
    glProgramUniform3fv(Program, GroundUniform_LightRadiance, 1, Frame.LightRadiance);
    glProgramUniform1f(Program, GroundUniform_ZNear, Frame.ZNear);
    glProgramUniform1f(Program, GroundUniform_ClearDepth, Frame.ClearDepth);
    glProgramUniform1f(Program, GroundUniform_DepthScale, Impl->DepthScale);
    glProgramUniform1f(Program, GroundUniform_TargetHeight, static_cast<float>(Impl->Height));

    glUseProgram(Program);
    glBindVertexArray(Impl->GroundVertexArray.Get());
    glDrawArrays(GL_TRIANGLES, 0, 3);
}

void FGpuRenderer::EndFrame()
{
    MOO_CHECK(Impl->bInFrame);
    Impl->bInFrame = false;
    Impl->bGroundAllowed = false;
    glBindVertexArray(0);
    glUseProgram(0);
    glDisable(GL_DEPTH_CLAMP);
    glDisable(GL_PROGRAM_POINT_SIZE);
    glDisable(GL_SCISSOR_TEST);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    Impl->Timer.End();
}

void FGpuRenderer::ReadbackRGB8(FImageRGB8& OutImage) const
{
    MOO_CHECK(!Impl->bInFrame);
    const int Width = Impl->Width;
    const int Height = Impl->Height;
    const size_t RowBytes = static_cast<size_t>(Width) * 3;

    std::vector<uint8_t> BottomUp(RowBytes * Height);
    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
    glPixelStorei(GL_PACK_ROW_LENGTH, 0);
    glPixelStorei(GL_PACK_SKIP_ROWS, 0);
    glPixelStorei(GL_PACK_SKIP_PIXELS, 0);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glGetTextureImage(Impl->ColorTarget.Get(), 0, GL_RGB, GL_UNSIGNED_BYTE, static_cast<GLsizei>(BottomUp.size()), BottomUp.data());

    //GL纹理行0为底部，输出约定行0为顶部
    OutImage.Width = Width;
    OutImage.Height = Height;
    OutImage.Pixels.resize(BottomUp.size());
    for(int Row = 0; Row < Height; Row++)
    {
        std::memcpy(&OutImage.Pixels[static_cast<size_t>(Row) * RowBytes], &BottomUp[static_cast<size_t>(Height - 1 - Row) * RowBytes], RowBytes);
    }
}

void FGpuRenderer::PresentToBackbuffer(int BackbufferWidth, int BackbufferHeight) const
{
    MOO_CHECK(!Impl->bInFrame);
    //最小化时帧缓冲尺寸为0，跳过呈现
    if(BackbufferWidth <= 0 || BackbufferHeight <= 0)
    {
        return;
    }
    glBlitNamedFramebuffer(Impl->Framebuffer.Get(), 0,
        0, 0, Impl->Width, Impl->Height,
        0, 0, BackbufferWidth, BackbufferHeight,
        GL_COLOR_BUFFER_BIT, GL_NEAREST);
}

FViewportTexture FGpuRenderer::GetViewportTexture() const
{
    return {Impl->ColorTarget.Get(), Impl->Width, Impl->Height, 0.f, 1.f, 1.f, 0.f};
}

double FGpuRenderer::GetLastGpuMs() const
{
    return Impl->Timer.GetLastMs();
}
