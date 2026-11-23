#ifndef MOOLAB_CPU_IMAGE_PRESENTER
#define MOOLAB_CPU_IMAGE_PRESENTER

// CPU参考路径的呈现器：把CPU帧缓冲（RGB float 0~255，行0为顶部）上传为纹理并全屏绘制到默认帧缓冲。
// 与FGpuRenderer无依赖，--renderer cpu 只需要窗口与本呈现器。实现位于 CpuImagePresenter.cpp（pimpl）。

#include <memory>
#include "ViewportTexture.h"

class IGlContextOwner;
struct FImageRGBFView;

class FCpuImagePresenter
{
public:
    // 需在Window之后构造、之前析构（GL资源依赖窗口的GL上下文）
    FCpuImagePresenter(const IGlContextOwner& Context, int Width, int Height);
    ~FCpuImagePresenter();

    FCpuImagePresenter(const FCpuImagePresenter&) = delete;
    FCpuImagePresenter& operator=(const FCpuImagePresenter&) = delete;

    // Image尺寸必须与构造时一致（Width*Height*3个float）
    void Present(const FImageRGBFView& Image, int BackbufferWidth, int BackbufferHeight);
    void UpdateViewport(const FImageRGBFView& Image);
    FViewportTexture GetViewportTexture() const;

private:
    struct FImpl;
    std::unique_ptr<FImpl> Impl;
};

#endif
