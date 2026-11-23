#include "GlRhi.h"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <vector>

#include "Core.h"
#include "GlContext.h"

void CheckGlContextReady(const IGlContextOwner& Context)
{
    MOO_CHECK_ALWAYS(Context.IsGlContextCurrent(), "renderer requires its owner's current GL context on this thread");
    MOO_CHECK_ALWAYS(glad_glGetIntegerv != nullptr && GLAD_GL_VERSION_4_6, "OpenGL 4.6 functions not loaded");
    GLint Major = 0, Minor = 0;
    glGetIntegerv(GL_MAJOR_VERSION, &Major);
    glGetIntegerv(GL_MINOR_VERSION, &Minor);
    MOO_CHECK_ALWAYS(Major > 4 || (Major == 4 && Minor >= 6), "OpenGL 4.6 required, got %d.%d", Major, Minor);
}

void GlDeleteBuffer(GLuint Id) { glDeleteBuffers(1, &Id); }
void GlDeleteVertexArray(GLuint Id) { glDeleteVertexArrays(1, &Id); }
void GlDeleteTexture(GLuint Id) { glDeleteTextures(1, &Id); }
void GlDeleteFramebuffer(GLuint Id) { glDeleteFramebuffers(1, &Id); }
void GlDeleteProgram(GLuint Id) { glDeleteProgram(Id); }
void GlDeleteQuery(GLuint Id) { glDeleteQueries(1, &Id); }

FGlBuffer CreateGlBuffer(GLsizeiptr Size, const void* Data, GLbitfield Flags)
{
    GLuint Id = 0;
    glCreateBuffers(1, &Id);
    glNamedBufferStorage(Id, Size, Data, Flags);
    return FGlBuffer(Id);
}

FGlVertexArray CreateGlVertexArray()
{
    GLuint Id = 0;
    glCreateVertexArrays(1, &Id);
    return FGlVertexArray(Id);
}

FGlTexture2D CreateGlTexture2D(GLenum InternalFormat, int Width, int Height)
{
    GLuint Id = 0;
    glCreateTextures(GL_TEXTURE_2D, 1, &Id);
    glTextureStorage2D(Id, 1, InternalFormat, Width, Height);
    glTextureParameteri(Id, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTextureParameteri(Id, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTextureParameteri(Id, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTextureParameteri(Id, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    return FGlTexture2D(Id);
}

FGlFramebuffer CreateGlFramebuffer()
{
    GLuint Id = 0;
    glCreateFramebuffers(1, &Id);
    return FGlFramebuffer(Id);
}

namespace
{
std::string ReadTextFile(const std::string& Path)
{
    std::ifstream File(Path, std::ios::binary);
    MOO_CHECK_ALWAYS(File.good(), "failed to open shader file '%s'", Path.c_str());
    std::ostringstream Stream;
    Stream << File.rdbuf();
    return Stream.str();
}

GLuint CompileShaderSource(GLenum Stage, const std::string& Source, const std::string& Description)
{
    const char* SourcePtr = Source.c_str();

    const GLuint Shader = glCreateShader(Stage);
    glShaderSource(Shader, 1, &SourcePtr, nullptr);
    glCompileShader(Shader);

    GLint bCompiled = GL_FALSE;
    glGetShaderiv(Shader, GL_COMPILE_STATUS, &bCompiled);
    if(bCompiled != GL_TRUE)
    {
        GLint LogLength = 0;
        glGetShaderiv(Shader, GL_INFO_LOG_LENGTH, &LogLength);
        std::vector<char> Log(static_cast<size_t>(LogLength > 1 ? LogLength : 1), '\0');
        glGetShaderInfoLog(Shader, static_cast<GLsizei>(Log.size()), nullptr, Log.data());
        MOO_CHECK_ALWAYS(false, "failed to compile shader %s:\n%s", Description.c_str(), Log.data());
    }
    return Shader;
}

GLuint CompileShader(GLenum Stage, const std::string& Path)
{
    return CompileShaderSource(Stage, ReadTextFile(Path), "'" + Path + "'");
}

// 主文件首行须为#version；公共源码拼在其后，再接主文件余文。用#line把源串号设为0=主文件、1=公共文件，
// 诊断中的“串号(行号)”可对应回各自文件行。不支持递归include。
GLuint CompileShaderWithCommon(GLenum Stage, const std::string& Path, const std::string& CommonPath)
{
    const std::string Main = ReadTextFile(Path);
    const std::string Common = ReadTextFile(CommonPath);
    const size_t FirstLineEnd = Main.find('\n');
    MOO_CHECK_ALWAYS(Main.compare(0, 8, "#version") == 0 && FirstLineEnd != std::string::npos,
                     "shader '%s' must start with a #version line to accept common source", Path.c_str());
    MOO_CHECK_ALWAYS(Common.find("#version") == std::string::npos, "common shader source '%s' must not contain #version", CommonPath.c_str());
    std::string Source = Main.substr(0, FirstLineEnd + 1);
    Source += "#line 1 1\n";
    Source += Common;
    Source += "\n#line 2 0\n";
    Source += Main.substr(FirstLineEnd + 1);
    return CompileShaderSource(Stage, Source, "'" + Path + "' (source 0) + common '" + CommonPath + "' (source 1)");
}
}

FGlProgram CreateGlProgramFromFiles(const std::string& VertexPath, const std::string& FragmentPath, const std::string& CommonFragmentPath)
{
    const GLuint VertexShader = CompileShader(GL_VERTEX_SHADER, VertexPath);
    const GLuint FragmentShader = CommonFragmentPath.empty() ? CompileShader(GL_FRAGMENT_SHADER, FragmentPath)
        : CompileShaderWithCommon(GL_FRAGMENT_SHADER, FragmentPath, CommonFragmentPath);

    FGlProgram Program(glCreateProgram());
    glAttachShader(Program.Get(), VertexShader);
    glAttachShader(Program.Get(), FragmentShader);
    glLinkProgram(Program.Get());
    //链接完成后着色器对象可删除，程序对象保留已链接的二进制
    glDetachShader(Program.Get(), VertexShader);
    glDetachShader(Program.Get(), FragmentShader);
    glDeleteShader(VertexShader);
    glDeleteShader(FragmentShader);

    GLint bLinked = GL_FALSE;
    glGetProgramiv(Program.Get(), GL_LINK_STATUS, &bLinked);
    if(bLinked != GL_TRUE)
    {
        GLint LogLength = 0;
        glGetProgramiv(Program.Get(), GL_INFO_LOG_LENGTH, &LogLength);
        std::vector<char> Log(static_cast<size_t>(LogLength > 1 ? LogLength : 1), '\0');
        glGetProgramInfoLog(Program.Get(), static_cast<GLsizei>(Log.size()), nullptr, Log.data());
        MOO_CHECK_ALWAYS(false, "failed to link program '%s' + '%s':\n%s", VertexPath.c_str(), FragmentPath.c_str(), Log.data());
    }
    return Program;
}

std::string GetGlShaderPath(const char* FileName)
{
#ifdef MOOLAB_SHADER_DIR
    return std::string(MOOLAB_SHADER_DIR) + "/" + FileName;
#else
    return std::string("Shaders/") + FileName;
#endif
}

FGlGpuTimer::FGlGpuTimer()
{
    for(int i = 0; i < RingSize; i++)
    {
        GLuint Ids[2] = {};
        glCreateQueries(GL_TIMESTAMP, 2, Ids);
        BeginQueries[i] = FGlQuery(Ids[0]);
        EndQueries[i] = FGlQuery(Ids[1]);
    }
}

void FGlGpuTimer::Begin()
{
    CollectAvailable();
    // 环满时跳过采样，未完成查询既不覆盖也不等待，计时不成为GPU同步点。
    bSampling = !bPending[WriteIndex];
    if(bSampling)
    {
        glQueryCounter(BeginQueries[WriteIndex].Get(), GL_TIMESTAMP);
    }
}

void FGlGpuTimer::End()
{
    if(!bSampling)
    {
        return;
    }
    bSampling = false;
    glQueryCounter(EndQueries[WriteIndex].Get(), GL_TIMESTAMP);
    bPending[WriteIndex] = true;
    WriteIndex = (WriteIndex + 1) % RingSize;
}

void FGlGpuTimer::CollectAvailable()
{
    //从最旧的一组开始按提交顺序读取，遇到未就绪的即停止（其后的更新结果也不会就绪）
    for(int Offset = 0; Offset < RingSize; Offset++)
    {
        const int Index = (WriteIndex + Offset) % RingSize;
        if(!bPending[Index])
        {
            continue;
        }
        GLint bAvailable = GL_FALSE;
        glGetQueryObjectiv(EndQueries[Index].Get(), GL_QUERY_RESULT_AVAILABLE, &bAvailable);
        if(bAvailable != GL_TRUE)
        {
            break;
        }
        GLuint64 BeginNs = 0, EndNs = 0;
        glGetQueryObjectui64v(BeginQueries[Index].Get(), GL_QUERY_RESULT, &BeginNs);
        glGetQueryObjectui64v(EndQueries[Index].Get(), GL_QUERY_RESULT, &EndNs);
        LastMs = static_cast<double>(EndNs - BeginNs) * 1e-6;
        bPending[Index] = false;
    }
}
