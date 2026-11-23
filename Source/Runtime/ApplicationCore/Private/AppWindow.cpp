#include "AppWindow.h"

#include <cstdio>
#include <cstdlib>
#include <filesystem>

#include <glad/glad.h>
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

#include "Core.h"

namespace
{
void OnGlfwError(int Code, const char* Description)
{
    std::fprintf(stderr, "[GLFW] error 0x%X: %s\n", Code, Description);
}

const char* GetGlDebugTypeName(GLenum Type)
{
    switch(Type)
    {
        case GL_DEBUG_TYPE_ERROR:               return "ERROR";
        case GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR: return "DEPRECATED";
        case GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR:  return "UNDEFINED";
        case GL_DEBUG_TYPE_PORTABILITY:         return "PORTABILITY";
        case GL_DEBUG_TYPE_PERFORMANCE:         return "PERFORMANCE";
        default:                                return "OTHER";
    }
}

const char* GetGlDebugSeverityName(GLenum Severity)
{
    switch(Severity)
    {
        case GL_DEBUG_SEVERITY_HIGH:   return "HIGH";
        case GL_DEBUG_SEVERITY_MEDIUM: return "MEDIUM";
        case GL_DEBUG_SEVERITY_LOW:    return "LOW";
        default:                       return "NOTIFICATION";
    }
}

void APIENTRY OnGlDebugMessage(GLenum Source, GLenum Type, GLuint Id, GLenum Severity,
                               GLsizei Length, const GLchar* Message, const void* UserParam)
{
    (void)Source;
    (void)Length;
    (void)UserParam;
    //通知级消息（如缓冲放置位置提示）量大且无害，不输出
    if(Severity == GL_DEBUG_SEVERITY_NOTIFICATION)
    {
        return;
    }
    std::fprintf(stderr, "[GL] type=%s severity=%s id=%u: %s\n",
                 GetGlDebugTypeName(Type), GetGlDebugSeverityName(Severity), Id, Message);
}
}

FAppWindow::FAppWindow(const FAppWindowDesc& Desc)
{
    glfwSetErrorCallback(OnGlfwError);
    MOO_CHECK_ALWAYS(glfwInit() == GLFW_TRUE, "glfwInit failed");

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GLFW_TRUE);
    glfwWindowHint(GLFW_VISIBLE, Desc.bHidden ? GLFW_FALSE : GLFW_TRUE);
#ifndef NDEBUG
    glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GLFW_TRUE);
#endif

    Window = glfwCreateWindow(Desc.Width, Desc.Height, Desc.Title, nullptr, nullptr);
    MOO_CHECK_ALWAYS(Window != nullptr, "failed to create OpenGL 4.6 core window (driver lacks GL 4.6?)");
    glfwMakeContextCurrent(Window);
    //粘滞键：同一次事件轮询内按下并抬起的短促按键至少被查询到一次，快捷键与ESC不会丢失
    glfwSetInputMode(Window, GLFW_STICKY_KEYS, GLFW_TRUE);
    //帧节拍由FFrameLimiter控制，关闭垂直同步避免两者叠加
    glfwSwapInterval(0);

    MOO_CHECK_ALWAYS(gladLoadGLLoader(reinterpret_cast<GLADloadproc>(glfwGetProcAddress)) != 0, "failed to load OpenGL functions");
    MOO_CHECK_ALWAYS(GLAD_GL_VERSION_4_6 != 0, "OpenGL 4.6 is required, got %d.%d", GLVersion.major, GLVersion.minor);

    GLint ContextFlags = 0;
    glGetIntegerv(GL_CONTEXT_FLAGS, &ContextFlags);
    if((ContextFlags & GL_CONTEXT_FLAG_DEBUG_BIT) != 0)
    {
        glEnable(GL_DEBUG_OUTPUT);
        //同步输出，保证消息紧随出错调用，便于定位
        glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
        glDebugMessageCallback(OnGlDebugMessage, nullptr);
    }

    std::fprintf(stdout, "[MooLab] OpenGL %s | %s\n",
                 reinterpret_cast<const char*>(glGetString(GL_VERSION)),
                 reinterpret_cast<const char*>(glGetString(GL_RENDERER)));

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    //不在资源目录生成imgui.ini
    auto& Io = ImGui::GetIO();
    Io.IniFilename = nullptr;
    Io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
    // 仅解析已安装字体；失败时Editor明确使用英文标签。
    if(const char* WindowsDirectory = std::getenv("WINDIR"))
    {
        for(const char* Name : {"msyh.ttc", "simhei.ttf", "simsun.ttc"})
        {
            const auto Path = std::filesystem::path(WindowsDirectory) / "Fonts" / Name;
            std::error_code Error;
            if(!std::filesystem::is_regular_file(Path, Error) || Error) continue;
            if(Io.Fonts->AddFontFromFileTTF(Path.string().c_str(), 18.f, nullptr, Io.Fonts->GetGlyphRangesChineseFull()))
            {
                bChineseFont = true;
                break;
            }
        }
    }
    ImGui::StyleColorsDark();
    MOO_CHECK_ALWAYS(ImGui_ImplGlfw_InitForOpenGL(Window, true), "ImGui GLFW backend initialization failed");
    MOO_CHECK_ALWAYS(ImGui_ImplOpenGL3_Init("#version 460 core"), "ImGui OpenGL3 backend initialization failed");
}

FAppWindow::~FAppWindow()
{
    //ImplOpenGL3会删除自身GL对象，须在上下文销毁前关闭；顺序与初始化相反
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(Window);
    glfwTerminate();
}

bool FAppWindow::IsGlContextCurrent() const
{
    return Window != nullptr && glfwGetCurrentContext() == Window;
}

void FAppWindow::BeginFrame()
{
    glfwPollEvents();
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
}

void FAppWindow::EndFrame()
{
    ImGui::Render();
    int Width = 0, Height = 0;
    glfwGetFramebufferSize(Window, &Width, &Height);
    if(Width > 0 && Height > 0)
    {
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glViewport(0, 0, Width, Height);
        glDisable(GL_SCISSOR_TEST);
        glDisable(GL_FRAMEBUFFER_SRGB);
        glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
        glClearColor(0.06f, 0.065f, 0.075f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(Window);
    }
}

bool FAppWindow::ShouldClose() const
{
    //直接查询按键状态，不安装自定义按键回调，避免与ImplGlfw的回调链交织；控件获得焦点时ESC同样退出
    return glfwWindowShouldClose(Window) != 0 || glfwGetKey(Window, GLFW_KEY_ESCAPE) == GLFW_PRESS;
}

void FAppWindow::GetFramebufferSize(int& OutWidth, int& OutHeight) const
{
    glfwGetFramebufferSize(Window, &OutWidth, &OutHeight);
}
