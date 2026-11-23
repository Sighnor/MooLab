#include "Input.h"
#include "AppWindow.h"
#include <cmath>
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

FRawInput PollWindowInput(const FAppWindow& AppWindow)
{
    FRawInput Raw;
    GLFWwindow* Window = AppWindow.Window;
    // 窗口启用粘滞键：失焦时也查询一遍以清除残留状态，再整体清空，重新聚焦不会补发旧按键。
    auto Down = [&](int Key) { return glfwGetKey(Window, Key) == GLFW_PRESS; };
    Raw.bW = Down(GLFW_KEY_W); Raw.bA = Down(GLFW_KEY_A);
    Raw.bS = Down(GLFW_KEY_S); Raw.bD = Down(GLFW_KEY_D);
    Raw.bLeft = Down(GLFW_KEY_LEFT); Raw.bRight = Down(GLFW_KEY_RIGHT);
    Raw.bUp = Down(GLFW_KEY_UP); Raw.bDown = Down(GLFW_KEY_DOWN);
    Raw.bSpace = Down(GLFW_KEY_SPACE); Raw.bPeriod = Down(GLFW_KEY_PERIOD);
    Raw.bFocused = glfwGetWindowAttrib(Window, GLFW_FOCUSED) == GLFW_TRUE;
    if(!Raw.bFocused) return FRawInput{};
    glfwGetCursorPos(Window, &Raw.MouseX, &Raw.MouseY);
    Raw.bMouseLeft = glfwGetMouseButton(Window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
    Raw.bMouseRight = glfwGetMouseButton(Window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;
    return Raw;
}
namespace
{
FGamepadStick Filter(FGamepadStick Stick)
{
    const float Magnitude = std::sqrt(Stick.X * Stick.X + Stick.Y * Stick.Y);
    if(!std::isfinite(Magnitude) || Magnitude <= 0.2f) return {};
    const float Response = Magnitude > 1.f ? 1.f : Magnitude * Magnitude;
    return {Stick.X / Magnitude * Response, Stick.Y / Magnitude * Response};
}
}
FRoutedInput FInputRouter::Route(const FRawInput& Raw, const FGamepadState& Gamepad, bool bAllowScene)
{
    FRoutedInput Result;
    if(bAllowScene && Raw.bFocused)
    {
        // 键盘和手柄原始贡献相加后才执行旧径向响应。
        Result.Move = Filter({Gamepad.LeftStick.X + float(Raw.bD) - float(Raw.bA),
                              Gamepad.LeftStick.Y + float(Raw.bS) - float(Raw.bW)});
        Result.Look = Filter({Gamepad.RightStick.X + float(Raw.bRight) - float(Raw.bLeft),
                              Gamepad.RightStick.Y + float(Raw.bDown) - float(Raw.bUp)});
        Result.bTogglePause = Raw.bSpace && !bPreviousSpace;
        Result.bSingleStep = Raw.bPeriod && !bPreviousPeriod;
    }
    // 门控关闭时也记录按键状态：在控件中按住的键切回视口后不算新按下。
    // 窗口失焦时原始输入整体清零，重新聚焦时仍按住的键会作为一次新按下。
    bPreviousSpace = Raw.bSpace;
    bPreviousPeriod = Raw.bPeriod;
    return Result;
}
