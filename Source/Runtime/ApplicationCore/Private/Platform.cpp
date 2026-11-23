#include "Platform.h"
#include "Core.h"

#include <algorithm>
#include <cmath>
#include <thread>

#if defined(_WIN32)
#include <windows.h>
#include <timeapi.h>
#include <Xinput.h>
// 旧版SDK头文件可能未定义该标志（Win10 1803引入）
#ifndef CREATE_WAITABLE_TIMER_HIGH_RESOLUTION
#define CREATE_WAITABLE_TIMER_HIGH_RESOLUTION 0x00000002
#endif
#endif

namespace
{
    // 摇杆原始值归一化到[-1, 1]；|v|<=0.1归零，与原raylib GetGamepadAxisMovement的防漂移阈值一致
    float NormalizeThumb(short Value)
    {
        const float Normalized = std::clamp(static_cast<float>(Value) / 32767.f, -1.f, 1.f);
        return std::fabs(Normalized) > 0.1f ? Normalized : 0.f;
    }
}

FGamepadState PollGamepad(int PlayerIndex)
{
    FGamepadState State;
#if defined(_WIN32)
    if(PlayerIndex < 0 || PlayerIndex >= XUSER_MAX_COUNT)
    {
        return State;
    }

    // 对未连接槽位调用XInputGetState开销可达毫秒级，断开后每秒只重试一次
    using Clock = std::chrono::steady_clock;
    static bool bLastConnected[XUSER_MAX_COUNT] = {};
    static Clock::time_point NextRetry[XUSER_MAX_COUNT] = {};
    const Clock::time_point Now = Clock::now();
    if(!bLastConnected[PlayerIndex] && Now < NextRetry[PlayerIndex])
    {
        return State;
    }

    XINPUT_STATE RawState = {};
    if(XInputGetState(static_cast<DWORD>(PlayerIndex), &RawState) != ERROR_SUCCESS)
    {
        bLastConnected[PlayerIndex] = false;
        NextRetry[PlayerIndex] = Now + std::chrono::seconds(1);
        return State;
    }
    bLastConnected[PlayerIndex] = true;

    // XInput的Y轴向上为正，取反以保持原raylib(GLFW)的向下为正约定
    State.bConnected = true;
    State.LeftStick.X = NormalizeThumb(RawState.Gamepad.sThumbLX);
    State.LeftStick.Y = -NormalizeThumb(RawState.Gamepad.sThumbLY);
    State.RightStick.X = NormalizeThumb(RawState.Gamepad.sThumbRX);
    State.RightStick.Y = -NormalizeThumb(RawState.Gamepad.sThumbRY);
#else
    (void)PlayerIndex;
#endif
    return State;
}

FFrameLimiter::FFrameLimiter(double TargetFps)
    : FramePeriod(std::chrono::steady_clock::duration::zero())
    , FrameStart(std::chrono::steady_clock::now())
{
    MOO_CHECK_ALWAYS(TargetFps > 0.0 && std::isfinite(TargetFps), "invalid target fps %f", TargetFps);
    FramePeriod = std::chrono::duration_cast<std::chrono::steady_clock::duration>(std::chrono::duration<double>(1.0 / TargetFps));
#if defined(_WIN32)
    // 高精度可等待计时器（Win10 1803+），精度约0.5ms，且不像timeBeginPeriod那样修改全局时钟分辨率
    WaitTimer = CreateWaitableTimerExW(nullptr, nullptr, CREATE_WAITABLE_TIMER_HIGH_RESOLUTION, TIMER_MODIFY_STATE | SYNCHRONIZE);
    if(WaitTimer == nullptr)
    {
        // 不支持时退回sleep_until，需将计时器精度提到1ms，否则按默认~15.6ms节拍量化
        bRaisedTimerResolution = timeBeginPeriod(1) == TIMERR_NOERROR;
        std::fprintf(stderr, "[MooLab] high-resolution waitable timer unavailable, fallback to sleep (timeBeginPeriod(1) %s)\n",
            bRaisedTimerResolution ? "ok" : "failed");
    }
#endif
}

FFrameLimiter::~FFrameLimiter()
{
#if defined(_WIN32)
    if(WaitTimer != nullptr)
    {
        CloseHandle(static_cast<HANDLE>(WaitTimer));
    }
    if(bRaisedTimerResolution)
    {
        timeEndPeriod(1);
    }
#endif
}

void FFrameLimiter::BeginFrame()
{
    FrameStart = std::chrono::steady_clock::now();
}

double FFrameLimiter::GetElapsedMs() const
{
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - FrameStart).count();
}

void FFrameLimiter::EndFrame()
{
    const std::chrono::steady_clock::time_point Deadline = FrameStart + FramePeriod;
    const std::chrono::steady_clock::duration Remaining = Deadline - std::chrono::steady_clock::now();
    if(Remaining <= std::chrono::steady_clock::duration::zero())
    {
        return;
    }
#if defined(_WIN32)
    if(WaitTimer != nullptr)
    {
        // 相对时间以100ns为单位、取负值
        LARGE_INTEGER DueTime;
        DueTime.QuadPart = -static_cast<LONGLONG>(std::chrono::duration_cast<std::chrono::nanoseconds>(Remaining).count() / 100);
        if(SetWaitableTimer(static_cast<HANDLE>(WaitTimer), &DueTime, 0, nullptr, nullptr, FALSE))
        {
            WaitForSingleObject(static_cast<HANDLE>(WaitTimer), INFINITE);
            return;
        }
    }
#endif
    std::this_thread::sleep_until(Deadline);
}
