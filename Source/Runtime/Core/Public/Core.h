#ifndef MOOLAB_CORE
#define MOOLAB_CORE

// Core层：断言与致命错误。本头文件只含宏与inline函数，可被任意编译单元包含（无ODR风险）

#include <cstdarg>
#include <cstdio>
#include <cstdlib>

// __restrict__ 为GCC/Clang写法，MSVC为__restrict
#if defined(_MSC_VER)
#define MOO_RESTRICT __restrict
#else
#define MOO_RESTRICT __restrict__
#endif

// printf格式串的编译期检查
#if defined(_MSC_VER)
#include <sal.h>
#define MOO_PRINTF_FORMAT_STRING _Printf_format_string_
#define MOO_PRINTF_ATTRIBUTE(FormatIndex, FirstArgIndex)
#else
#define MOO_PRINTF_FORMAT_STRING
#define MOO_PRINTF_ATTRIBUTE(FormatIndex, FirstArgIndex) __attribute__((format(printf, FormatIndex, FirstArgIndex)))
#endif

[[noreturn]] inline void MooFatalError(const char* File, int Line, const char* Expr, MOO_PRINTF_FORMAT_STRING const char* Format, ...) MOO_PRINTF_ATTRIBUTE(4, 5);
[[noreturn]] inline void MooFatalError(const char* File, int Line, const char* Expr, MOO_PRINTF_FORMAT_STRING const char* Format, ...)
{
    std::fprintf(stderr, "[MooLab Fatal] %s(%d): check '%s' failed: ", File, Line, Expr);
    va_list Args;
    va_start(Args, Format);
    std::vfprintf(stderr, Format, Args);
    va_end(Args);
    std::fprintf(stderr, "\n");
    std::fflush(stderr);
    std::abort();
}

// MOO_CHECK：开发期不变式，Release(NDEBUG)下裁掉，语义同 UE check
#ifdef NDEBUG
#define MOO_CHECK(Expr) ((void)0)
#else
#define MOO_CHECK(Expr) do { if(!(Expr)) { MooFatalError(__FILE__, __LINE__, #Expr, "invariant violated"); } } while(0)
#endif

// MOO_CHECK_ALWAYS(Expr, Format, ...)：所有配置均生效，用于资源缺失等外部输入错误，
// Release下也必须给出明确报错而非空指针崩溃；消息为printf格式
#define MOO_CHECK_ALWAYS(Expr, ...) do { if(!(Expr)) { MooFatalError(__FILE__, __LINE__, #Expr, __VA_ARGS__); } } while(0)

#endif
