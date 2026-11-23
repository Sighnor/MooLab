#ifndef MOOLAB_GL_RHI
#define MOOLAB_GL_RHI

// OpenGL 4.6 后端的RAII薄封装（DSA）。注意：这里只是GL对象的所有权管理，并非可替换后端的RHI抽象层。
// 暴露GLuint等GL类型，只允许Renderer的GL实现文件（GpuRenderer.cpp、CpuImagePresenter.cpp）及Tools/Verify下的GL探针包含，Sandbox/Editor/Launch不得包含。
// 所有对象只能在GL上下文有效期间创建与析构（见AppWindow.h生命周期约束）。

#include <glad/glad.h>

#include <string>

class IGlContextOwner;
// 必须在创建任何GL对象（包括计时查询）之前验证调用线程的当前上下文。
void CheckGlContextReady(const IGlContextOwner& Context);

using FGlDeleteFn = void (*)(GLuint);

// 仅可移动的GL对象句柄，析构时调用DeleteFn释放
template<FGlDeleteFn DeleteFn>
class TGlHandle
{
public:
    TGlHandle() = default;
    explicit TGlHandle(GLuint InId) : Id(InId) {}
    ~TGlHandle() { Reset(); }

    TGlHandle(const TGlHandle&) = delete;
    TGlHandle& operator=(const TGlHandle&) = delete;

    TGlHandle(TGlHandle&& Other) noexcept : Id(Other.Id) { Other.Id = 0; }
    TGlHandle& operator=(TGlHandle&& Other) noexcept
    {
        if(this != &Other)
        {
            Reset();
            Id = Other.Id;
            Other.Id = 0;
        }
        return *this;
    }

    void Reset()
    {
        if(Id != 0)
        {
            DeleteFn(Id);
            Id = 0;
        }
    }

    GLuint Get() const { return Id; }

private:
    GLuint Id = 0;
};

void GlDeleteBuffer(GLuint Id);
void GlDeleteVertexArray(GLuint Id);
void GlDeleteTexture(GLuint Id);
void GlDeleteFramebuffer(GLuint Id);
void GlDeleteProgram(GLuint Id);
void GlDeleteQuery(GLuint Id);

using FGlBuffer = TGlHandle<&GlDeleteBuffer>;
using FGlVertexArray = TGlHandle<&GlDeleteVertexArray>;
using FGlTexture2D = TGlHandle<&GlDeleteTexture>;
using FGlFramebuffer = TGlHandle<&GlDeleteFramebuffer>;
using FGlProgram = TGlHandle<&GlDeleteProgram>;
using FGlQuery = TGlHandle<&GlDeleteQuery>;

// 不可变存储缓冲；Flags为glNamedBufferStorage标志（如GL_DYNAMIC_STORAGE_BIT）
FGlBuffer CreateGlBuffer(GLsizeiptr Size, const void* Data, GLbitfield Flags);
FGlVertexArray CreateGlVertexArray();
// 单层mip的不可变纹理，最近邻、clamp到边缘
FGlTexture2D CreateGlTexture2D(GLenum InternalFormat, int Width, int Height);
FGlFramebuffer CreateGlFramebuffer();
// 读取并编译链接顶点/片元着色器文件，失败时打印info log并终止。
// CommonFragmentPath非空时，把该文件（不含#version）拼在片元主文件首行#version之后、其余正文之前；不是通用include。
FGlProgram CreateGlProgramFromFiles(const std::string& VertexPath, const std::string& FragmentPath,
                                    const std::string& CommonFragmentPath = std::string());
// 着色器文件的完整路径：编译期宏MOOLAB_SHADER_DIR（源码Shaders目录）+ 文件名
std::string GetGlShaderPath(const char* FileName);

// GPU计时：每帧在命令流首尾写GL_TIMESTAMP，环形使用多组查询避免读回阻塞；
// 只读取已可用的最新结果，因此GetLastMs()通常滞后若干帧，结果未就绪前返回0
class FGlGpuTimer
{
public:
    FGlGpuTimer();

    void Begin();
    void End();
    double GetLastMs() const { return LastMs; }

private:
    static constexpr int RingSize = 4;
    FGlQuery BeginQueries[RingSize];
    FGlQuery EndQueries[RingSize];
    bool bPending[RingSize] = {};
    bool bSampling = false;
    int WriteIndex = 0;
    double LastMs = 0.0;

    void CollectAvailable();
};

#endif
