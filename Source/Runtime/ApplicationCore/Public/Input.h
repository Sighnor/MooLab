#pragma once
#include "Platform.h"

class FAppWindow;
struct FRawInput
{
    bool bFocused = false;
    bool bW = false, bA = false, bS = false, bD = false;
    bool bLeft = false, bRight = false, bUp = false, bDown = false;
    bool bSpace = false, bPeriod = false;
    double MouseX = 0., MouseY = 0.;
    bool bMouseLeft = false, bMouseRight = false;
};
struct FRoutedInput
{
    FGamepadStick Move, Look;
    bool bTogglePause = false;
    bool bSingleStep = false;
};
FRawInput PollWindowInput(const FAppWindow& Window);
class FInputRouter
{
public:
    FRoutedInput Route(const FRawInput& Raw, const FGamepadState& Gamepad, bool bAllowScene);
private:
    bool bPreviousSpace = false;
    bool bPreviousPeriod = false;
};
