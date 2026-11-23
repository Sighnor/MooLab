#include "Application.h"
#include "Core.h"
#include "AppWindow.h"
#include "Platform.h"
#include "SandboxScene.h"
#include "GpuRenderer.h"
#include "CpuImagePresenter.h"
#include "EditorUi.h"
#include "Input.h"
#include <algorithm>
#include <climits>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>

namespace
{
struct FAppOptions
{
    int MaxFrames = -1;
    bool bUnlimitedFps = false;
    bool bHidden = false;
    bool bGpuRenderer = true;
    ESceneMode Mode = ESceneMode::Physics;
    std::vector<std::pair<int, std::string>> Captures;
};
FAppOptions ParseAppOptions(int ArgCount, char** Arguments)
{
    auto ParseCount = [](const char* Text, int& Value)
    {
        char* End = nullptr;
        const long Parsed = std::strtol(Text, &End, 10);
        if(End == Text || *End != '\0' || Parsed < 0 || Parsed > INT_MAX) return false;
        Value = static_cast<int>(Parsed);
        return true;
    };
    FAppOptions Options;
    for(int Index = 1; Index < ArgCount; ++Index)
    {
        const char* Option = Arguments[Index];
        if(std::strcmp(Option, "--frames") == 0)
        {
            MOO_CHECK_ALWAYS(Index + 1 < ArgCount, "incomplete option: --frames");
            const char* Value = Arguments[++Index];
            MOO_CHECK_ALWAYS(ParseCount(Value, Options.MaxFrames), "--frames expects a non-negative integer, got '%s'", Value);
        }
        else if(std::strcmp(Option, "--capture-frame") == 0)
        {
            MOO_CHECK_ALWAYS(Index + 2 < ArgCount, "incomplete option: --capture-frame");
            int Frame = 0;
            const char* Value = Arguments[++Index];
            MOO_CHECK_ALWAYS(ParseCount(Value, Frame), "--capture-frame expects a non-negative frame index, got '%s'", Value);
            std::error_code Error;
            const auto Path = std::filesystem::absolute(Arguments[++Index], Error);
            MOO_CHECK_ALWAYS(!Error, "invalid capture path: %s", Error.message().c_str());
            Options.Captures.emplace_back(Frame, Path.string());
        }
        else if(std::strcmp(Option, "--renderer") == 0 || std::strcmp(Option, "--scene") == 0)
        {
            MOO_CHECK_ALWAYS(Index + 1 < ArgCount, "incomplete option: %s", Option);
            const char* Value = Arguments[++Index];
            if(std::strcmp(Option, "--renderer") == 0)
            {
                MOO_CHECK_ALWAYS(std::strcmp(Value, "gpu") == 0 || std::strcmp(Value, "cpu") == 0,
                                 "invalid renderer '%s': expected gpu or cpu", Value);
                Options.bGpuRenderer = std::strcmp(Value, "gpu") == 0;
            }
            else
            {
                MOO_CHECK_ALWAYS(std::strcmp(Value, "physics") == 0 || std::strcmp(Value, "motion-matching") == 0,
                                 "invalid scene '%s': expected physics or motion-matching", Value);
                Options.Mode = std::strcmp(Value, "physics") == 0 ? ESceneMode::Physics : ESceneMode::MotionMatching;
            }
        }
        else if(std::strcmp(Option, "--hidden") == 0) Options.bHidden = true;
        else if(std::strcmp(Option, "--unlimited-fps") == 0) Options.bUnlimitedFps = true;
        else std::cerr << "Unknown option: " << Option << '\n';
    }
    return Options;
}
}

int RunApplication(int ArgCount, char** Arguments)
{
    const FAppOptions Options = ParseAppOptions(ArgCount, Arguments);
    FFrameLimiter Limiter(60.);
    FAppWindow Window({1600, 900, "MooLab Editor", Options.bHidden});
    FSandboxScene Scene;
    FSceneCommand InitialCommand;
    InitialCommand.bSetMode = true;
    InitialCommand.Mode = Options.Mode;
    if(!Scene.ApplyCommand(InitialCommand)) { std::cerr << Scene.GetSnapshot().Message << '\n'; return 1; }
    std::unique_ptr<FGpuRenderer> Gpu;
    std::unique_ptr<FCpuImagePresenter> Cpu;
    if(Options.bGpuRenderer) { Gpu = std::make_unique<FGpuRenderer>(Window, 720, 720); Scene.InitializeGpu(*Gpu); }
    else Cpu = std::make_unique<FCpuImagePresenter>(Window, 720, 720);
    FEditorUi Editor(Window.HasChineseFont());
    FInputRouter Router;
    using FClock = std::chrono::steady_clock;
    const auto Start = FClock::now();
    auto Previous = Start;
    FFrameStats Stats;
    double MaxInterval = 0.;
    bool bCaptureFailed = false;
    while(!Window.ShouldClose() && (Options.MaxFrames < 0 || Stats.AppFrame < static_cast<uint64_t>(Options.MaxFrames)))
    {
        Limiter.BeginFrame();
        const auto Now = FClock::now();
        if(Stats.AppFrame > 0)
        {
            Stats.IntervalMs = std::chrono::duration<double, std::milli>(Now - Previous).count();
            MaxInterval = std::max(MaxInterval, Stats.IntervalMs);
        }
        Previous = Now;
        Window.BeginFrame();
        const auto Texture = Gpu ? Gpu->GetViewportTexture() : Cpu->GetViewportTexture();
        auto Ui = Editor.Build(Scene.GetSettings(), Scene.GetSnapshot(), Stats, Texture);
        const auto Raw = Options.bHidden ? FRawInput{} : PollWindowInput(Window);
        const auto Gamepad = Options.bHidden ? FGamepadState{} : PollGamepad(0);
        const auto Routed = Router.Route(Raw, Gamepad, !Options.bHidden && Ui.bAllowSceneInput);
        Ui.Command.bTogglePause |= Routed.bTogglePause;
        Ui.Command.bSingleStep |= Routed.bSingleStep;
        if(Ui.bSettingsChanged) Scene.SetSettings(Ui.Settings);
        Scene.ApplyCommand(Ui.Command);
        Scene.Tick({Routed.Move.X, Routed.Move.Y, Routed.Look.X, Routed.Look.Y});
        if(Gpu) { Scene.BuildRenderData(*Gpu); Stats.GpuMs = Gpu->GetLastGpuMs(); }
        else Cpu->UpdateViewport(Scene.RenderCpu());
        const bool bCapture = std::any_of(Options.Captures.begin(), Options.Captures.end(),
            [&](const auto& Capture) { return static_cast<uint64_t>(Capture.first) == Stats.AppFrame; });
        if(bCapture)
        {
            FImageRGB8 Image;
            if(Gpu) Gpu->ReadbackRGB8(Image); else Scene.ReadCpuRGB8(Image);
            cv::Mat Rgb(Image.Height, Image.Width, CV_8UC3, Image.Pixels.data()), Bgr;
            cv::cvtColor(Rgb, Bgr, cv::COLOR_RGB2BGR);
            for(const auto& Capture : Options.Captures)
            {
                if(static_cast<uint64_t>(Capture.first) != Stats.AppFrame) continue;
                try
                {
                    if(!cv::imwrite(Capture.second, Bgr)) { bCaptureFailed = true; std::cerr << "Failed to write capture: " << Capture.second << '\n'; }
                }
                catch(const cv::Exception& Error) { bCaptureFailed = true; std::cerr << "Failed to write capture: " << Capture.second << ": " << Error.what() << '\n'; }
            }
        }
        Window.EndFrame();
        Stats.CostMs = Limiter.GetElapsedMs();
        std::ostringstream Log;
        Log << " Time: " << Stats.AppFrame * 0.0166667f << " Frame: " << Stats.AppFrame << " Cost: " << Stats.CostMs
            << "ms Interval: " << Stats.IntervalMs << "ms Gpu: " << Stats.GpuMs << "ms\n";
        std::cout << Log.str();
        ++Stats.AppFrame;
        if(!Options.bUnlimitedFps) Limiter.EndFrame();
    }
    if(Stats.AppFrame > 0)
        std::cout << "Summary: " << Stats.AppFrame << " frames, avg frame interval "
                  << std::chrono::duration<double, std::milli>(FClock::now() - Start).count() / Stats.AppFrame
                  << "ms, max frame interval " << MaxInterval << "ms\n";
    for(const auto& Capture : Options.Captures)
        if(static_cast<uint64_t>(Capture.first) >= Stats.AppFrame)
        {
            bCaptureFailed = true;
            std::cerr << "Warning: capture frame " << Capture.first << " not reached (stopped at frame " << Stats.AppFrame << "): " << Capture.second << '\n';
        }
    return bCaptureFailed ? 1 : 0;
}
