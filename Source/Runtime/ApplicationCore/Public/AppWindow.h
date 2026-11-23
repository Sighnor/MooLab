#ifndef MOOLAB_APP_WINDOW
#define MOOLAB_APP_WINDOW

// 应用窗口层：GLFW窗口 + OpenGL 4.6 core上下文 + Dear ImGui生命周期，替代原OpenCV highgui窗口与trackbar。
// 实现位于 AppWindow.cpp；本头文件不暴露glad/GLFW（仅前置声明），可被Launch等非GL模块安全包含。
//
// 生命周期约束：所有GL资源（FGpuRenderer、FCpuImagePresenter）必须在窗口之后构造、之前析构，
// 二者仅依赖IGlContextOwner契约，构造时校验本窗口的上下文是否为调用线程的当前上下文。
// 每帧顺序：BeginFrame（处理消息、开始ImGui帧）→ ImGui控件与场景绘制 → 呈现到默认帧缓冲 → EndFrame（绘制ImGui、交换缓冲）。

#include "GlContext.h"
#include "Input.h"

struct GLFWwindow;

struct FAppWindowDesc
{
    int Width = 720;
    int Height = 720;
    const char* Title = "MOOLAB";
    // 隐藏窗口（无人值守验证），渲染与截图走离屏帧缓冲，不受影响
    bool bHidden = false;
};

class FAppWindow : public IGlContextOwner
{
public:
    explicit FAppWindow(const FAppWindowDesc& Desc);
    ~FAppWindow();

    FAppWindow(const FAppWindow&) = delete;
    FAppWindow& operator=(const FAppWindow&) = delete;

    bool IsGlContextCurrent() const override;
    void BeginFrame();
    // 绘制ImGui到默认帧缓冲并交换缓冲
    void EndFrame();
    // 点击关闭按钮或按下ESC
    bool ShouldClose() const;
    // 默认帧缓冲的像素尺寸，最小化时为0
    void GetFramebufferSize(int& OutWidth, int& OutHeight) const;

public:
    bool HasChineseFont() const { return bChineseFont; }
private:
    friend FRawInput PollWindowInput(const FAppWindow& Window);
    GLFWwindow* Window = nullptr;
    bool bChineseFont = false;
};

#endif
