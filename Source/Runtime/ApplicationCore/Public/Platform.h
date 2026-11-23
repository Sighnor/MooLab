#ifndef MOOLAB_PLATFORM
#define MOOLAB_PLATFORM

// 平台层：手柄输入与帧节拍，替代原先仅用于读手柄的raylib。
// 实现位于 Platform.cpp（独立编译单元，隔离 <windows.h> 的宏污染）。
// 本头文件只放声明，不得包含带非inline定义的存量头文件（vec.hpp等），否则多编译单元下ODR重复定义。

#include <chrono>

// 摇杆值范围[-1, 1]；Y轴向下为正，与原raylib(GLFW)手柄约定一致。
// 仅做与原raylib相同的每轴|v|<=0.1归零（防漂移），径向死区与响应曲线由调用方处理
struct FGamepadStick
{
    float X = 0.f;
    float Y = 0.f;
};

struct FGamepadState
{
    bool bConnected = false;
    FGamepadStick LeftStick;
    FGamepadStick RightStick;
};

// 读取指定玩家(0~3)的手柄状态，未连接时返回全零。内部有重试退避状态，只允许在单一（主）线程调用
FGamepadState PollGamepad(int PlayerIndex);

// 固定帧率节拍器：BeginFrame记录帧起点，EndFrame等待到下一帧时刻。
// 优先使用高精度可等待计时器；系统不支持时在生命周期内将进程计时器精度提到1ms（timeBeginPeriod），
// 与原raylib InitWindow的行为一致
class FFrameLimiter
{
public:
    explicit FFrameLimiter(double TargetFps);
    ~FFrameLimiter();

    FFrameLimiter(const FFrameLimiter&) = delete;
    FFrameLimiter& operator=(const FFrameLimiter&) = delete;

    void BeginFrame();
    // 自BeginFrame起经过的墙钟时间（毫秒）
    double GetElapsedMs() const;
    // 等待至本帧结束时刻（BeginFrame + 1/TargetFps），超时则立即返回
    void EndFrame();

private:
    std::chrono::steady_clock::duration FramePeriod;
    std::chrono::steady_clock::time_point FrameStart;
    void* WaitTimer = nullptr;
    bool bRaisedTimerResolution = false;
};

#endif
