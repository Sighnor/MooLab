#ifndef MOOLAB_GPU_RENDERER
#define MOOLAB_GPU_RENDERER

// OpenGL 4.6 GPU渲染路径：顶点变换与逐像素着色在GPU执行，着色逻辑逐行移植自CPU参考实现（shader.hpp）。
// 实现位于 GpuRenderer.cpp（pimpl，头文件不暴露GL类型）；数据接口与坐标约定见 RenderTypes.h。
//
// 每帧调用顺序：BeginFrame → [DrawGround] → DrawMesh/DrawPixelSquares（按CPU绘制顺序）→ EndFrame → [ReadbackRGB8] → PresentToBackbuffer。
// 后续GPU蒙皮/FK等compute阶段应插在BeginFrame之前（写顶点缓冲后加glMemoryBarrier）；
// 网格各属性使用独立缓冲（不交错），便于compute直接写位置/法线。

#include <memory>

#include "RenderTypes.h"

class FGpuRenderer
{
public:
    // 需在Window之后构造、之前析构（GL资源依赖窗口的GL上下文）
    FGpuRenderer(const IGlContextOwner& Context, int Width, int Height);
    ~FGpuRenderer();

    FGpuRenderer(const FGpuRenderer&) = delete;
    FGpuRenderer& operator=(const FGpuRenderer&) = delete;

    FMeshId CreateMesh(const FMeshDataView& Mesh);
    // 更新已创建网格的全部位置/法线，顶点数必须不变；UV/索引保持静态，不在BeginFrame/EndFrame之间调用。
    void UpdateMeshVertices(FMeshId Mesh, const float* Positions, const float* Normals, int VertexCount);
    // 按CPU存储原样上传为RGB32F（值0~255），着色器取值与CPU逐位一致
    FTextureId CreateTextureRGBF(const FImageRGBFView& Image);

    void BeginFrame(const FFrameParams& Params);
    // 程序化地面（逐像素解析求交，写入view depth供后续网格深度测试）。须为帧内第一个绘制且至多一次，
    // 帧外、重复或在DrawMesh/DrawPixelSquares之后调用视为契约错误；参数不可渲染时调用方应直接跳过。
    void DrawGround(const FGroundParams& Ground);
    // Model为16个float列主序；Pbr着色需提供PbrParams，Light着色可传nullptr。
    // 像素映射与覆盖范围复刻CPU光栅（含整体上移1行、最后一列与图像底行不着色），截图可与CPU参考逐像素直接对拍
    void DrawMesh(FMeshId Mesh, const float Model[16], EShading Shading, const FPbrParams* PbrParams);
    // 在像素(x, y)（行0为顶部）处绘制边长2*HalfSize+1的实心方块，不做深度测试；ColorRGB取值0~1
    void DrawPixelSquares(const int* PixelXY, int Count, int HalfSize, const float ColorRGB[3]);
    void EndFrame();

    // 读回本帧离屏颜色（不含ImGui），行0为顶部
    void ReadbackRGB8(FImageRGB8& OutImage) const;
    // 将离屏颜色拉伸blit到默认帧缓冲
    void PresentToBackbuffer(int BackbufferWidth, int BackbufferHeight) const;
    // 最近一次可用的GPU帧耗时（毫秒），结果滞后若干帧，未就绪时为0
    double GetLastGpuMs() const;
    FViewportTexture GetViewportTexture() const;

private:
    struct FImpl;
    std::unique_ptr<FImpl> Impl;
};

#endif
