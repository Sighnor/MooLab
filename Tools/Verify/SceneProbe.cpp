#include "SandboxScene.h"
#include "AppWindow.h"
#include "GpuRenderer.h"
#include "PbrShading.h"
#include "Input.h"
#include "FootContacts.h"
#include "TwoBoneIk.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>

namespace
{
constexpr float Pi = 3.14159265358979f;
constexpr float StepSeconds = 0.0166667f;
constexpr int Seams[] = {0, 281, 562, 13153, 25744, 39622, 53500};

void Require(bool Condition, const char* Name)
{
    if(!Condition) { std::cerr << "FAIL: " << Name << '\n'; throw Name; }
    std::cout << "PASS: " << Name << '\n';
}
void Require(bool Condition, const std::string& Name) { Require(Condition, Name.c_str()); }
int BudgetFailures = 0;
void Budget(bool Condition, const char* Name, const std::string& Detail)
{
    // 继续收集后续证据；退出时用码3区分预算失败与硬检查失败(码1)。
    if(!Condition) ++BudgetFailures;
    std::cout << (Condition ? "BUDGET PASS: " : "BUDGET FAIL: ") << Name << " (" << Detail << ")\n";
}
void Budget(bool Condition, const std::string& Name, const std::string& Detail) { Budget(Condition, Name.c_str(), Detail); }
void Save(const FImageRGB8& Image, const std::string& Path)
{
    cv::Mat Rgb(Image.Height, Image.Width, CV_8UC3, const_cast<uint8_t*>(Image.Pixels.data())), Bgr;
    cv::cvtColor(Rgb, Bgr, cv::COLOR_RGB2BGR);
    Require(cv::imwrite(Path, Bgr), "capture written");
}
float Wrap(float Angle) { return std::remainder(Angle, 2.f * Pi); }
float Degrees(float Radians) { return Radians * 180.f / Pi; }
int SegmentOf(int Frame)
{
    for(int Index = 1; Index < 7; ++Index) if(Frame < Seams[Index]) return Index - 1;
    return -1;
}
bool InClipSegment(int Frame, int Start)
{
    return Start >= 1 && Start + 75 < 53500 && Frame >= Start - 1 && Frame <= Start + 75
        && SegmentOf(Start - 1) == SegmentOf(Start + 75);
}
bool SnapshotFinite(const FSceneSnapshot& Snapshot)
{
    bool bFinite = std::isfinite(Snapshot.RootYaw) && std::isfinite(Snapshot.RootSpeed);
    for(float Value : Snapshot.RootPosition) bFinite &= std::isfinite(Value);
    for(float Value : Snapshot.CameraPosition) bFinite &= std::isfinite(Value);
    for(const auto& Bone : Snapshot.Bones)
    {
        for(float Value : Bone.Position) bFinite &= std::isfinite(Value);
        for(float Value : Bone.Rotation) bFinite &= std::isfinite(Value);
    }
    return bFinite;
}
// CameraYaw为弧度；相机跟随Hips，不能用根骨减相机位置估计视线。
std::array<float, 2> CameraForward(const FSceneSnapshot& Snapshot)
{
    return {std::sin(Snapshot.CameraYaw), std::cos(Snapshot.CameraYaw)};
}
// 10步水平平均速度，抑制步态周期内的根摆动。
struct FMotionTracker
{
    std::vector<std::array<float, 2>> History;
    std::vector<float> Speeds;
    void Push(const FSceneSnapshot& Snapshot)
    {
        History.push_back({Snapshot.RootPosition[0], Snapshot.RootPosition[2]});
        Speeds.push_back(Snapshot.RootSpeed);
    }
    float MeanSpeed() const
    {
        if(Speeds.size() < 10) return std::numeric_limits<float>::infinity();
        float Sum = 0.f;
        for(size_t Index = Speeds.size() - 10; Index < Speeds.size(); ++Index) Sum += Speeds[Index];
        return Sum / 10.f;
    }
    std::array<float, 2> Velocity() const
    {
        if(History.size() < 11) return {0.f, 0.f};
        const auto& A = History[History.size() - 11];
        const auto& B = History.back();
        return {(B[0] - A[0]) / (10.f * StepSeconds), (B[1] - A[1]) / (10.f * StepSeconds)};
    }
};
float Speed(const std::array<float, 2>& Velocity) { return std::sqrt(Velocity[0] * Velocity[0] + Velocity[1] * Velocity[1]); }
float AngleBetween(const std::array<float, 2>& A, const std::array<float, 2>& B)
{
    return std::fabs(std::atan2(A[0] * B[1] - A[1] * B[0], A[0] * B[0] + A[1] * B[1]));
}
FSceneInput MoveInput(const std::array<float, 2>& World, const std::array<float, 2>& Forward)
{
    // 世界方向 -> 相机局部输入：Right = cross(Forward, up) = (-Fz, Fx)。
    const float Ahead = World[0] * Forward[0] + World[1] * Forward[1];
    const float Side = World[0] * -Forward[1] + World[1] * Forward[0];
    return {Side, -Ahead, 0.f, 0.f};
}
void Tick(FSandboxScene& Scene, FMotionTracker& Tracker, const FSceneInput& Input)
{
    Scene.Tick(Input);
    Tracker.Push(Scene.GetSnapshot());
}
void EnterMatching(FSandboxScene& Scene)
{
    FSceneCommand Command;
    Command.bSetMode = true;
    Command.Mode = ESceneMode::MotionMatching;
    Require(Scene.ApplyCommand(Command), "enter motion matching");
}

void ProbeRouter()
{
    FInputRouter Router;
    FRawInput Raw;
    Raw.bFocused = true;
    Raw.bW = Raw.bS = true;
    auto Routed = Router.Route(Raw, {}, true);
    Require(Routed.Move.X == 0.f && Routed.Move.Y == 0.f, "W and S cancel");
    Raw.bW = Raw.bS = false;
    Raw.bA = Raw.bD = true;
    Routed = Router.Route(Raw, {}, true);
    Require(Routed.Move.X == 0.f && Routed.Move.Y == 0.f, "A and D cancel");
    Raw = {};
    Raw.bFocused = true;
    Raw.bW = Raw.bD = true;
    Routed = Router.Route(Raw, {}, true);
    Require(std::fabs(std::hypot(Routed.Move.X, Routed.Move.Y) - 1.f) < 1e-6f && Routed.Move.Y < 0.f && Routed.Move.X > 0.f,
            "diagonal keys saturate");
    FGamepadState Pad;
    Pad.LeftStick = {0.15f, 0.f};
    Raw = {};
    Raw.bFocused = true;
    Routed = Router.Route(Raw, Pad, true);
    Require(Routed.Move.X == 0.f, "deadzone after sum");
    Raw.bD = true;
    Pad.LeftStick = {-0.9f, 0.f};
    Routed = Router.Route(Raw, Pad, true);
    Require(Routed.Move.X == 0.f, "keyboard and pad summed before filter");
    Raw = {};
    Raw.bFocused = true;
    Raw.bW = true;
    Raw.bSpace = true;
    Routed = Router.Route(Raw, {}, false);
    Require(Routed.Move.Y == 0.f && !Routed.bTogglePause, "blocked input routes nothing");
    Routed = Router.Route(Raw, {}, true);
    Require(!Routed.bTogglePause, "held space consumed while blocked");
    Raw.bSpace = false;
    Router.Route(Raw, {}, true);
    Raw.bSpace = true;
    Require(Router.Route(Raw, {}, true).bTogglePause && !Router.Route(Raw, {}, true).bTogglePause, "space edge triggers once");
    Raw.bFocused = false;
    Require(Router.Route(Raw, {}, true).Move.Y == 0.f, "unfocused input cleared");
}

void ProbeMatching(FSandboxScene& Scene)
{
    auto Snapshot = Scene.GetSnapshot();
    std::cout << "MM candidates=" << Snapshot.CandidateCount << " excluded=" << Snapshot.ExcludedCount << '\n';
    Require(Snapshot.bMatchingAvailable && Snapshot.CandidateCount > 0, "motion matching available");
    EnterMatching(Scene);
    Snapshot = Scene.GetSnapshot();
    Require(Snapshot.Mode == ESceneMode::MotionMatching && Snapshot.NextFrame == 0 && SnapshotFinite(Snapshot)
            && Snapshot.SourceFrame >= 0, "matching reset state");

    // 暂停零推进、单步+1。
    FSceneCommand Pause;
    Pause.bTogglePause = true;
    Scene.ApplyCommand(Pause);
    for(int Step = 0; Step < 10; ++Step) Scene.Tick({0.f, -1.f, 0.f, 0.f});
    const auto Paused = Scene.GetSnapshot();
    Require(Paused.SourceFrame == Snapshot.SourceFrame && Paused.RootPosition == Snapshot.RootPosition, "matching pause freezes");
    FSceneCommand Step;
    Step.bSingleStep = true;
    Scene.ApplyCommand(Step);
    Scene.Tick({0.f, -1.f, 0.f, 0.f});
    Scene.Tick({0.f, -1.f, 0.f, 0.f});
    Require(Scene.GetSnapshot().NextFrame == 1, "matching single step once");
    Scene.ApplyCommand(Pause);

    // 零输入、非有限输入。
    Scene.Reset();
    for(int Index = 0; Index < 240; ++Index) Scene.Tick({});
    Snapshot = Scene.GetSnapshot();
    Require(SnapshotFinite(Snapshot), "zero input finite");
    std::cout << "zero input speed after 240 steps=" << Snapshot.RootSpeed << '\n';
    const float NaN = std::numeric_limits<float>::quiet_NaN();
    for(int Index = 0; Index < 30; ++Index) Scene.Tick({NaN, NaN, 0.f, 0.f});
    Require(SnapshotFinite(Scene.GetSnapshot()), "nonfinite input ignored");

    // 起步、停步。
    Scene.Reset();
    FMotionTracker Tracker;
    const auto Forward = CameraForward(Scene.GetSnapshot());
    int StartStep = -1;
    float StartError = 0.f;
    for(int Index = 1; Index <= 120; ++Index)
    {
        Tick(Scene, Tracker, MoveInput(Forward, Forward));
        const auto Velocity = Tracker.Velocity();
        if(StartStep < 0 && Speed(Velocity) >= 1.f && AngleBetween(Velocity, Forward) <= Pi / 6.f)
        { StartStep = Index; StartError = Degrees(AngleBetween(Velocity, Forward)); }
    }
    Budget(StartStep > 0 && StartStep <= 60, "W start <=60 steps speed>=1 error<=30deg",
           "step=" + std::to_string(StartStep) + " error=" + std::to_string(StartError)
           + " speed@120=" + std::to_string(Speed(Tracker.Velocity())));
    std::vector<float> StopSpeeds;
    int FirstBelow = -1;
    float MinFirst60 = std::numeric_limits<float>::infinity();
    float SpeedAt60 = 0.f;
    for(int Index = 1; Index <= 180; ++Index)
    {
        Tick(Scene, Tracker, {});
        const float Current = Tracker.MeanSpeed();
        StopSpeeds.push_back(Current);
        if(Index <= 60) MinFirst60 = std::min(MinFirst60, Current);
        if(Index == 60) SpeedAt60 = Current;
        if(FirstBelow < 0 && Current <= 0.2f) FirstBelow = Index;
        if(Index % 10 == 0)
        {
            const auto D = Scene.GetSnapshot();
            std::cout << "  stop " << Index << " meanSpeed=" << Current << " src=" << D.SourceFrame << " clip=" << D.ClipStart
                      << " n=" << D.SearchCount << " cost=" << D.SearchCost << '\n';
        }
    }
    // 在前60步内寻找任一满足点，并检查其后完整60步，而非只检查首次短暂降速。
    int StopStep = -1;
    for(int Index = 0; Index < 60; ++Index)
    {
        bool bHolds = true;
        for(int Offset = 0; Offset <= 60; ++Offset) bHolds &= StopSpeeds[Index + Offset] <= 0.2f;
        if(bHolds) { StopStep = Index + 1; break; }
    }
    Budget(StopStep > 0, "release stops <=60 steps and holds 60",
           "step=" + std::to_string(StopStep) + " firstBelow=" + std::to_string(FirstBelow)
           + " minFirst60=" + std::to_string(MinFirst60) + " speed@60=" + std::to_string(SpeedAt60));

    // 90°与180°转向。
    const std::array<float, 2> Right = {-Forward[1], Forward[0]};
    const std::array<float, 2> Back = {-Forward[0], -Forward[1]};
    for(int Case = 0; Case < 2; ++Case)
    {
        Scene.Reset();
        FMotionTracker Turn;
        for(int Index = 0; Index < 120; ++Index) Tick(Scene, Turn, MoveInput(Forward, Forward));
        const auto Target = Case == 0 ? Right : Back;
        float LastYaw = Scene.GetSnapshot().RootYaw;
        float Turned = 0.f;
        float BestError = 180.f;
        float HeadingError = 180.f;
        int BackwardSteps = 0;
        int Reached = -1;
        for(int Index = 1; Index <= 90; ++Index)
        {
            Tick(Scene, Turn, MoveInput(Target, Forward));
            const auto Now = Scene.GetSnapshot();
            Turned += Wrap(Now.RootYaw - LastYaw);
            LastYaw = Now.RootYaw;
            const auto Velocity = Turn.Velocity();
            const std::array<float, 2> Facing = {std::sin(Now.RootYaw), std::cos(Now.RootYaw)};
            const float Error = Degrees(AngleBetween(Velocity, Target));
            HeadingError = Degrees(AngleBetween(Facing, Target));
            if(Speed(Velocity) > 0.5f && Velocity[0] * Facing[0] + Velocity[1] * Facing[1] < 0.f) ++BackwardSteps;
            if(Index % 10 == 0)
            {
                const auto V = Turn.Velocity();
                std::cout << "  turn" << Case << ' ' << Index << " yaw=" << Degrees(Now.RootYaw) << " velDir="
                          << Degrees(std::atan2(V[0], V[1])) << " v=" << Speed(V) << " src=" << Now.SourceFrame
                          << " clip=" << Now.ClipStart << " n=" << Now.SearchCount << '\n';
            }
            BestError = std::min(BestError, Error);
            if(Reached < 0 && (Case == 0 ? Error <= 30.f && Speed(Turn.Velocity()) > 0.3f : std::fabs(Degrees(Turned)) >= 150.f))
                Reached = Index;
        }
        if(Case == 0)
            Budget(Reached > 0, "90deg turn error<=30deg within 90 steps",
                   "step=" + std::to_string(Reached) + " best=" + std::to_string(BestError));
        else
            Budget(Reached > 0, "180deg turn >=150deg within 90 steps",
                   "step=" + std::to_string(Reached) + " turned=" + std::to_string(Degrees(Turned)));
        std::cout << "turn" << Case << " headingError@90=" << HeadingError
                  << " backwardSteps(speed>0.5)=" << BackwardSteps << '\n';
    }

    // 相机持续环绕时期望方向跨越±π。
    Scene.Reset();
    int Crossings = 0;
    float MaxWrappedDelta = 0.f;
    float MaxRawDelta = 0.f;
    float PreviousYaw = Scene.GetSnapshot().RootYaw;
    bool bFinite = true;
    for(int Index = 0; Index < 1200; ++Index)
    {
        Scene.Tick({0.f, -1.f, 0.05f, 0.f});
        const auto Now = Scene.GetSnapshot();
        bFinite &= SnapshotFinite(Now);
        Crossings += std::fabs(Now.RootYaw - PreviousYaw) > Pi;
        MaxRawDelta = std::max(MaxRawDelta, std::fabs(Degrees(Now.RootYaw - PreviousYaw)));
        MaxWrappedDelta = std::max(MaxWrappedDelta, std::fabs(Degrees(Wrap(Now.RootYaw - PreviousYaw))));
        PreviousYaw = Now.RootYaw;
    }
    std::cout << "yaw wrap crossings=" << Crossings << " maxWrappedDelta=" << MaxWrappedDelta << '\n';
    Require(bFinite && Crossings > 0, "yaw crosses +-pi finite");
    Require(MaxRawDelta > 20.f, "negative control unwrapped yaw rejected");
    Budget(MaxWrappedDelta <= 20.f, "yaw wrap step <=20deg", "max=" + std::to_string(MaxWrappedDelta));
    Require(InClipSegment(100, 100) && !InClipSegment(281, 281) && !InClipSegment(200, 100),
            "clip predicate rejects seam and out-of-window negative controls");
}

struct FLongRun
{
    std::vector<float> Trace;
    float MaxDisplacement = 0.f;
    float MaxYawDelta = 0.f;
    std::string MaxYawWhere;
    int Searches = 0;
    double SearchMs = 0.;
    double TickMs = 0.;
};
FLongRun RunRandom(FSandboxScene& Scene, int Steps, uint32_t Seed = 20240531u)
{
    FLongRun Run;
    Scene.Reset();
    auto Next = [&]() { Seed = Seed * 1664525u + 1013904223u; return Seed >> 8; };
    FSceneInput Input;
    int Remaining = 0;
    auto Previous = Scene.GetSnapshot();
    int LastSearch = Previous.SearchCount;
    const float NaN = std::numeric_limits<float>::quiet_NaN();
    const auto Begin = std::chrono::steady_clock::now();
    for(int Index = 0; Index < Steps; ++Index)
    {
        if(Remaining-- <= 0)
        {
            // 分段常量输入：八方向、停止、模拟量、偶发非有限值与相机环绕。
            Remaining = 10 + Next() % 150;
            const uint32_t Kind = Next() % 16;
            Input = {};
            if(Kind < 9)
            {
                static const float Dirs[9][2] = {{0, -1}, {1, -1}, {1, 0}, {1, 1}, {0, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {0, 0}};
                Input.MoveX = Dirs[Kind][0];
                Input.MoveY = Dirs[Kind][1];
            }
            else if(Kind < 14)
            {
                Input.MoveX = (Next() % 2001) / 1000.f - 1.f;
                Input.MoveY = (Next() % 2001) / 1000.f - 1.f;
            }
            else if(Kind == 14) Input.MoveX = Input.MoveY = NaN;
            if(Next() % 4 == 0) Input.LookX = (Next() % 2001) / 1000.f - 1.f;
        }
        Scene.Tick(Input);
        const auto Now = Scene.GetSnapshot();
        if(!SnapshotFinite(Now)) { std::cerr << "nonfinite at step " << Index << '\n'; Require(false, "random run finite"); }
        const bool bInWindow = InClipSegment(Now.SourceFrame, Now.ClipStart);
        if(!bInWindow)
        {
            std::cerr << "source " << Now.SourceFrame << " clip " << Now.ClipStart << " at step " << Index << '\n';
            Require(false, "source frame inside clip segment");
        }
        const float DX = Now.RootPosition[0] - Previous.RootPosition[0];
        const float DZ = Now.RootPosition[2] - Previous.RootPosition[2];
        Run.MaxDisplacement = std::max(Run.MaxDisplacement, std::sqrt(DX * DX + DZ * DZ));
        const float YawDelta = std::fabs(Degrees(Wrap(Now.RootYaw - Previous.RootYaw)));
        if(YawDelta > Run.MaxYawDelta)
        {
            Run.MaxYawDelta = YawDelta;
            Run.MaxYawWhere = "step=" + std::to_string(Index) + " src=" + std::to_string(Previous.SourceFrame) + "->"
                + std::to_string(Now.SourceFrame) + " clip=" + std::to_string(Now.ClipStart) + " blending="
                + std::to_string(Now.bBlending) + " searchStep=" + std::to_string(Now.SearchFrame);
        }
        if(Now.SearchCount != LastSearch) { ++Run.Searches; Run.SearchMs += Now.SearchMs; LastSearch = Now.SearchCount; }
        Run.Trace.insert(Run.Trace.end(), {Now.RootPosition[0], Now.RootPosition[1], Now.RootPosition[2], Now.RootYaw,
                                           static_cast<float>(Now.SourceFrame)});
        Previous = Now;
    }
    Run.TickMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - Begin).count() / Steps;
    return Run;
}
// 运动确定性摘要：按位FNV-1a累积root、源帧与全部胶囊位姿；只读快照，不触碰渲染。
struct FMotionDigest
{
    uint64_t Hash = 1469598103934665603ull;
    void Bytes(const void* Data, size_t Size)
    {
        const auto* Byte = static_cast<const unsigned char*>(Data);
        for(size_t Index = 0; Index < Size; ++Index) { Hash ^= Byte[Index]; Hash *= 1099511628211ull; }
    }
    void Push(const FSceneSnapshot& Snapshot)
    {
        Bytes(Snapshot.RootPosition.data(), sizeof(Snapshot.RootPosition));
        Bytes(&Snapshot.RootYaw, sizeof(Snapshot.RootYaw));
        Bytes(&Snapshot.SourceFrame, sizeof(Snapshot.SourceFrame));
        for(const auto& Bone : Snapshot.Bones)
        {
            Bytes(Bone.Position.data(), sizeof(Bone.Position));
            Bytes(Bone.Rotation.data(), sizeof(Bone.Rotation));
            Bytes(&Bone.bContact, sizeof(Bone.bContact));
        }
    }
};
uint64_t PhysicsDigest(FSandboxScene& Scene, int Steps)
{
    FMotionDigest Digest;
    Scene.Reset();
    for(int Index = 0; Index < Steps; ++Index) { Scene.Tick({}); Digest.Push(Scene.GetSnapshot()); }
    return Digest.Hash;
}
uint64_t MatchingDigest(FSandboxScene& Scene, int Steps)
{
    FMotionDigest Digest;
    Scene.Reset();
    uint32_t Seed = 20240531u;
    auto Next = [&]() { Seed = Seed * 1664525u + 1013904223u; return Seed >> 8; };
    FSceneInput Input;
    for(int Index = 0; Index < Steps; ++Index)
    {
        if(Index % 60 == 0)
        {
            static const float Dirs[9][2] = {{0, -1}, {1, -1}, {1, 0}, {1, 1}, {0, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {0, 0}};
            const uint32_t Kind = Next() % 9;
            Input = {Dirs[Kind][0], Dirs[Kind][1], (Next() % 3) * 0.5f - 0.5f, 0.f};
        }
        Scene.Tick(Input);
        Digest.Push(Scene.GetSnapshot());
    }
    return Digest.Hash;
}
void PrintDigests(FSandboxScene& Scene)
{
    FSceneCommand Mode;
    Mode.bSetMode = true;
    Mode.Mode = ESceneMode::Physics;
    Scene.ApplyCommand(Mode);
    const uint64_t Physics = PhysicsDigest(Scene, 600);
    Mode.Mode = ESceneMode::MotionMatching;
    Require(Scene.ApplyCommand(Mode), "digest enters matching");
    const uint64_t Matching = MatchingDigest(Scene, 3000);
    std::cout << "motion digest physics600=" << std::hex << Physics << " matching3000=" << Matching << std::dec << '\n';
    Mode.Mode = ESceneMode::Physics;
    Scene.ApplyCommand(Mode);
}

// ---- 程序化地面与独立材质 ----
bool Near(float A, double B, double Tolerance = 1e-5) { return std::fabs(static_cast<double>(A) - B) <= Tolerance; }
bool NearVector(const float A[3], double X, double Y, double Z, double Tolerance = 1e-5)
{
    return Near(A[0], X, Tolerance) && Near(A[1], Y, Tolerance) && Near(A[2], Z, Tolerance);
}
FFrameParams ProbeFrame(float X, float Y, float Z)
{
    FFrameParams Frame;
    for(int Index = 0; Index < 16; Index += 5) Frame.View[Index] = Frame.Proj[Index] = 1.f;
    Frame.CameraPos[0] = X;
    Frame.CameraPos[1] = Y;
    Frame.CameraPos[2] = Z;
    Frame.ZNear = 0.05f;
    Frame.ClearDepth = 100.f;
    return Frame;
}
// 把像素(0,0)的射线替换为指定方向，并用微小像素增量使盒滤波足迹取下限
FGroundParams Aimed(FGroundParams Ground, float X, float Y, float Z)
{
    Ground.RayTopLeft[0] = X;
    Ground.RayTopLeft[1] = Y;
    Ground.RayTopLeft[2] = Z;
    for(int Index = 0; Index < 3; ++Index) Ground.RayStepX[Index] = Ground.RayStepY[Index] = 0.f;
    Ground.RayStepX[0] = 1e-5f;
    Ground.RayStepY[2] = 1e-5f;
    return Ground;
}
bool Hit(const FGroundParams& Ground, const float Camera[3], FGroundSample& Sample)
{
    return EvaluateGroundPixel(Ground, Camera, 0.05f, 100.f, 0, 0, Sample);
}

void ProbeGroundGeometry()
{
    FFrameParams Frame = ProbeFrame(0.f, 1.f, 0.f);
    const double Horizontal[4] = {0., 1., 0., 0.};
    FGroundParams Ground;
    Require(BuildGroundParams(Frame, Horizontal, 0.5f, 720, 720, Ground), "ground params build for horizontal plane");
    Require(Near(Ground.RayTopLeft[0], -719. / 720., 1e-6) && Near(Ground.RayTopLeft[1], 717. / 720., 1e-6) && Ground.RayTopLeft[2] == -1.f,
            "top-left 720^2 pixel ndc (-719/720, 717/720)");
    Require(NearVector(Ground.RayStepX, 2. / 720., 0., 0., 1e-7) && NearVector(Ground.RayStepY, 0., -2. / 720., 0., 1e-7),
            "pixel ray increments 2/W, -2/H");
    Require(NearVector(Ground.AxisU, 1., 0., 0., 0.) && NearVector(Ground.AxisV, 0., 0., 1., 0.) && NearVector(Ground.Origin, 0., 0., 0., 0.),
            "horizontal basis U=+X V=+Z origin 0");

    // 非单位View独立锚点：camera-to-world = Y(+90度)*X(-30度)，观察射线(-7/16,1/6,-1)
    FFrameParams Rotated = ProbeFrame(0.f, 1.f, 0.f);
    const double Cos30 = std::sqrt(3.) / 2.;
    Rotated.View[0] = 0.f; Rotated.View[1] = -0.5f; Rotated.View[2] = static_cast<float>(Cos30);
    Rotated.View[4] = 0.f; Rotated.View[5] = static_cast<float>(Cos30); Rotated.View[6] = 0.5f;
    Rotated.View[8] = -1.f; Rotated.View[9] = 0.f; Rotated.View[10] = 0.f;
    Rotated.Proj[0] = 2.f;
    Rotated.Proj[5] = 3.f;
    FGroundParams RotatedGround;
    const double Qx = -1. / 12. - Cos30;
    const double Qy = Cos30 / 6. - 0.5;
    const double Qz = 7. / 16.;
    Require(BuildGroundParams(Rotated, Horizontal, 0.5f, 8, 6, RotatedGround)
            && NearVector(RotatedGround.RayTopLeft, Qx, Qy, Qz, 1e-7), "rotated camera ray uses View transpose with independent XYZ anchor");
    Require(NearVector(RotatedGround.RayStepX, 0., 0., -1. / 8., 1e-7)
            && NearVector(RotatedGround.RayStepY, 1. / 18., -std::sqrt(3.) / 18., 0., 1e-7), "rotated rectangular camera pixel-step independent anchors");
    FGroundSample RotatedSample;
    const double ExpectedDepth = -1. / Qy;
    Require(Hit(RotatedGround, Rotated.CameraPos, RotatedSample) && Near(RotatedSample.Depth, ExpectedDepth)
            && NearVector(RotatedSample.Position, ExpectedDepth * Qx, 0., ExpectedDepth * Qz), "rotated camera intersection position and view depth independent anchor");

    // 独立几何锚点：C=(0,1,0)、q=(0,-0.25,-1)
    FGroundSample Sample;
    Require(Hit(Aimed(Ground, 0.f, -0.25f, -1.f), Frame.CameraPos, Sample) && Near(Sample.Depth, 4.)
            && NearVector(Sample.Position, 0., 0., -4.), "anchor y=0 t=4 P=(0,0,-4)");
    const double Slanted[4] = {0., 1., 1., 0.};
    FGroundParams Slope;
    Require(BuildGroundParams(Frame, Slanted, 0.5f, 720, 720, Slope), "ground params build for slanted plane");
    Require(Hit(Aimed(Slope, 0.f, -0.25f, -1.f), Frame.CameraPos, Sample) && Near(Sample.Depth, 0.8)
            && NearVector(Sample.Position, 0., 0.8, -0.8), "anchor y+z=0 t=0.8 P=(0,0.8,-0.8)");
    const double Expected = 1. / std::sqrt(2.);
    Require(NearVector(Slope.Plane, 0., Expected, Expected, 1e-7) && Slope.Plane[3] == 0.f
            && NearVector(Slope.AxisU, 1., 0., 0., 1e-7) && NearVector(Slope.AxisV, 0., -Expected, Expected, 1e-7),
            "slanted plane normalized basis U=+X V=U x n");

    // 默认Physics平面 (0,1,0,0.05) 即 y=-0.05，不叠加ContactEpsilon
    const double Physics[4] = {0., 1., 0., 0.05};
    FGroundParams Contact;
    Require(BuildGroundParams(Frame, Physics, 0.5f, 720, 720, Contact) && NearVector(Contact.Origin, 0., -0.05, 0., 1e-7)
            && Hit(Aimed(Contact, 0.f, -1.f, 0.f), Frame.CameraPos, Sample) && Near(Sample.Position[1], -0.05) && Near(Sample.Depth, 1.05),
            "default physics plane sign y=-0.05");
    const double Scaled[4] = {0., 4., 0., 0.2};
    FGroundParams ScaledGround;
    Require(BuildGroundParams(Frame, Scaled, 0.5f, 720, 720, ScaledGround) && NearVector(ScaledGround.Origin, 0., -0.05, 0., 1e-7),
            "plane coefficients normalized by true length");

    // near拒绝、100接受、超远/背面/平行/背向不画
    const float Far[3] = {0.f, 1.5625f, 0.f};
    Require(Hit(Aimed(Ground, 0.f, -0.015625f, -1.f), Far, Sample) && Sample.Depth == 100.f, "t == ClearDepth(100) accepted");
    const float Beyond[3] = {0.f, 1.5626f, 0.f};
    Require(!Hit(Aimed(Ground, 0.f, -0.015625f, -1.f), Beyond, Sample), "t > ClearDepth rejected");
    const float AtNear[3] = {0.f, 0.05f, 0.f};
    const float PastNear[3] = {0.f, 0.0501f, 0.f};
    Require(!Hit(Aimed(Ground, 0.f, -1.f, 0.f), AtNear, Sample) && Hit(Aimed(Ground, 0.f, -1.f, 0.f), PastNear, Sample),
            "t == ZNear rejected, just beyond accepted");
    const float Below[3] = {0.f, -1.f, 0.f};
    FFrameParams BelowFrame = ProbeFrame(0.f, -1.f, 0.f);
    FGroundParams Unused;
    Require(!BuildGroundParams(BelowFrame, Horizontal, 0.5f, 720, 720, Unused) && !Hit(Aimed(Ground, 0.f, 1.f, 0.f), Below, Sample),
            "camera behind plane not drawn");
    Require(!Hit(Aimed(Ground, 0.f, 0.f, -1.f), Frame.CameraPos, Sample) && !Hit(Aimed(Ground, 0.f, 0.1f, -1.f), Frame.CameraPos, Sample),
            "parallel and upward rays rejected");
    const double Degenerate[4] = {0., 0., 0., 1.};
    const double NonFinite[4] = {0., std::numeric_limits<double>::quiet_NaN(), 0., 0.};
    Require(!BuildGroundParams(Frame, Degenerate, 0.5f, 720, 720, Unused) && !BuildGroundParams(Frame, NonFinite, 0.5f, 720, 720, Unused),
            "degenerate/non-finite plane skipped");

    // 世界固定格子：同一格中心沿U移动0.5m翻色、1m回原色；负坐标按floor奇偶
    auto CoverageAt = [&](const float Camera[3], float X, float Z)
    {
        FGroundSample Result;
        const bool bHit = Hit(Aimed(Ground, X - Camera[0], -Camera[1], Z - Camera[2]), Camera, Result);
        return bHit ? Result.Coverage : -1.f;
    };
    const float Camera[3] = {0.f, 1.f, 0.f};
    Require(CoverageAt(Camera, 0.25f, 0.25f) > 0.99f && CoverageAt(Camera, 0.75f, 0.25f) < 0.01f && CoverageAt(Camera, 1.25f, 0.25f) > 0.99f,
            "cell center flips at +0.5m and returns at +1m along U");
    Require(CoverageAt(Camera, -0.25f, 0.25f) < 0.01f && CoverageAt(Camera, -0.25f, -0.25f) > 0.99f, "negative coordinates parity");
    const float Moved[3] = {3.f, 2.f, -1.f};
    Require(CoverageAt(Moved, 0.25f, 0.25f) > 0.99f && CoverageAt(Moved, 0.75f, 0.25f) < 0.01f, "camera move keeps world anchor in same cell");
    FGroundParams MovedGround;
    BuildGroundParams(ProbeFrame(3.f, 2.f, -1.f), Horizontal, 0.5f, 720, 720, MovedGround);
    Require(std::equal(Ground.Plane, Ground.Plane + 4, MovedGround.Plane) && std::equal(Ground.Origin, Ground.Origin + 3, MovedGround.Origin)
            && std::equal(Ground.AxisU, Ground.AxisU + 3, MovedGround.AxisU) && std::equal(Ground.AxisV, Ground.AxisV + 3, MovedGround.AxisV),
            "grid frame depends only on plane");
    // 盒滤波中等足迹独立锚点：v=.5处保持+1，u在格边界取均值0，在u=.75,w=1时均值.5
    auto FilteredAt = [&](float X)
    {
        FGroundParams Filtered = Aimed(Ground, X, -1.f, 0.25f);
        Filtered.RayStepX[0] = 0.5f;
        Filtered.RayStepY[2] = 0.5f;
        FGroundSample Result;
        const bool bHit = Hit(Filtered, Camera, Result);
        return bHit ? Result.Coverage : -1.f;
    };
    Require(Near(FilteredAt(0.25f), 1., 1e-6) && Near(FilteredAt(0.5f), 0.5, 1e-6)
            && Near(FilteredAt(0.375f), 0.75, 1e-6) && Near(FilteredAt(-0.125f), 0.25, 1e-6), "box filter medium footprint independent anchors (including negative coordinate)");
    FGroundParams TwoAxes = Aimed(Ground, 0.375f, -1.f, 0.25f);
    for(int Index = 0; Index < 3; ++Index) TwoAxes.RayStepX[Index] = TwoAxes.RayStepY[Index] = 0.f;
    TwoAxes.RayStepX[0] = TwoAxes.RayStepY[0] = 0.25f;
    // v足迹取下限会产生float相消，误差仍受独立位置锚点1e-5约束
    FGroundSample TwoAxesSample;
    Require(Hit(TwoAxes, Camera, TwoAxesSample) && Near(TwoAxesSample.Coverage, 0.75, 1e-5), "footprint sums both pixel axes (not max or L2)");
    // 非零dt锚点：t=2、u=.9、v=.5；dP=(.775,0,-.125)，wu=1.55、wv=.25，覆盖率35/62
    FGroundParams Derivative = Aimed(Ground, 0.225f, -0.5f, 0.125f);
    Derivative.RayStepX[0] = 0.5f;
    Derivative.RayStepX[1] = -0.25f;
    for(float& Value : Derivative.RayStepY) Value = 0.f;
    FGroundSample DerivativeSample;
    Require(Hit(Derivative, Camera, DerivativeSample) && Near(DerivativeSample.Coverage, 35. / 62., 1e-6), "analytic dt/dP footprint independent anchor 35/62");
    // 有效中心 t=1，但半像素上边界为平行射线；未退化时格中心覆盖率应为1，故能区分该分支
    FGroundParams TouchHorizon = Aimed(Ground, 0.25f, -1.f, 0.25f);
    TouchHorizon.RayStepX[0] = 0.f;
    TouchHorizon.RayStepX[1] = 2.f;
    for(float& Value : TouchHorizon.RayStepY) Value = 0.f;
    FGroundSample HorizonSample;
    Require(Hit(TouchHorizon, Camera, HorizonSample) && HorizonSample.Depth == 1.f && HorizonSample.Coverage == 0.5f,
            "valid center footprint touching parallel singularity averages");
    TouchHorizon.RayStepX[1] = -2.f;
    Require(Hit(TouchHorizon, Camera, HorizonSample) && HorizonSample.Coverage == 0.5f, "horizon margin uses absolute pixel-step normal");
    TouchHorizon.RayStepX[1] = 1.5f;
    Require(Hit(TouchHorizon, Camera, HorizonSample) && HorizonSample.Coverage == 1.f, "half-pixel footprint below singularity stays filtered");
    TouchHorizon.RayStepX[1] = 0.f;
    TouchHorizon.RayStepY[1] = 2.f;
    Require(Hit(TouchHorizon, Camera, HorizonSample) && HorizonSample.Coverage == 0.5f, "parallel singularity also checks vertical pixel-step");
    // 远处足迹过宽退为平均色
    FGroundParams Wide = Aimed(Ground, 0.25f, -1.f, 0.25f);
    Wide.RayStepX[0] = Wide.RayStepY[2] = 100.f;
    FGroundSample WideSample;
    Require(Hit(Wide, Camera, WideSample) && WideSample.Coverage == 0.5f, "footprint >= 64 cells averages");
}

struct FDiffStats
{
    int Changed = 0;
    int Within3 = 0;
    int Over32 = 0;
    int Total = 0;
};
FDiffStats Compare(const FImageRGB8& A, const FImageRGB8& B)
{
    if(A.Width != B.Width || A.Height != B.Height) Require(false, "compared images have same size");
    FDiffStats Stats;
    Stats.Total = A.Width * A.Height;
    for(size_t Pixel = 0; Pixel < A.Pixels.size(); Pixel += 3)
    {
        int MaxDifference = 0;
        for(int Channel = 0; Channel < 3; ++Channel)
            MaxDifference = std::max(MaxDifference, std::abs(A.Pixels[Pixel + Channel] - B.Pixels[Pixel + Channel]));
        Stats.Changed += MaxDifference > 0;
        Stats.Within3 += MaxDifference <= 3;
        Stats.Over32 += MaxDifference > 32;
    }
    return Stats;
}
// 与MooImageDiff同一门槛：最大通道差<=3的像素>=99%，差>32的像素<=0.5%
bool Parity(const FImageRGB8& Gpu, const FImageRGB8& Cpu)
{
    const FDiffStats Stats = Compare(Gpu, Cpu);
    return Stats.Within3 >= 0.99 * Stats.Total && Stats.Over32 <= 0.005 * Stats.Total;
}
struct FFramePair
{
    FImageRGB8 Gpu;
    FImageRGB8 Cpu;
};
FFramePair Render(FSandboxScene& Scene, FGpuRenderer& Renderer)
{
    FFramePair Pair;
    Scene.BuildRenderData(Renderer);
    Renderer.ReadbackRGB8(Pair.Gpu);
    Scene.RenderCpu();
    Scene.ReadCpuRGB8(Pair.Cpu);
    return Pair;
}
void SetPaused(FSandboxScene& Scene, bool bPaused)
{
    if(Scene.GetSnapshot().bPaused == bPaused) return;
    FSceneCommand Toggle;
    Toggle.bTogglePause = true;
    Scene.ApplyCommand(Toggle);
}
void SetMode(FSandboxScene& Scene, ESceneMode Mode)
{
    FSceneCommand Command;
    Command.bSetMode = true;
    Command.Mode = Mode;
    Require(Scene.ApplyCommand(Command), Mode == ESceneMode::Physics ? "set physics mode" : "set matching mode");
}
bool SameMaterial(const FSurfaceMaterial& A, const FSurfaceMaterial& B) { return A.Metallic == B.Metallic && A.Roughness == B.Roughness; }
bool SameSurfaceSettings(const FSceneSettings& A, const FSceneSettings& B)
{
    bool bSame = SameMaterial(A.CharacterMaterial, B.CharacterMaterial) && SameMaterial(A.Ground.Material, B.Ground.Material)
        && A.Ground.bVisible == B.Ground.bVisible && A.Ground.CellSize == B.Ground.CellSize
        && A.Ground.MatchingHeight == B.Ground.MatchingHeight && A.ContactPlane == B.ContactPlane
        && A.bShowCapsules == B.bShowCapsules && A.bShowCharacter == B.bShowCharacter && A.bShowContacts == B.bShowContacts;
    for(int Index = 0; Index < SceneCapsuleCount; ++Index) bSame &= SameMaterial(A.CapsuleMaterials[Index], B.CapsuleMaterials[Index]);
    return bSame;
}

void ProbeSettingsValidation(FSandboxScene& Scene)
{
    using FField = std::function<float&(FSceneSettings&)>;
    const std::vector<std::pair<FField, std::pair<float, float>>> Fields = {
        {[](FSceneSettings& S) -> float& { return S.CharacterMaterial.Metallic; }, {0.f, 1.f}},
        {[](FSceneSettings& S) -> float& { return S.CharacterMaterial.Roughness; }, {0.05f, 1.f}},
        {[](FSceneSettings& S) -> float& { return S.CapsuleMaterials[7].Metallic; }, {0.f, 1.f}},
        {[](FSceneSettings& S) -> float& { return S.CapsuleMaterials[SceneCapsuleCount - 1].Roughness; }, {0.05f, 1.f}},
        {[](FSceneSettings& S) -> float& { return S.Ground.Material.Metallic; }, {0.f, 1.f}},
        {[](FSceneSettings& S) -> float& { return S.Ground.Material.Roughness; }, {0.05f, 1.f}},
        {[](FSceneSettings& S) -> float& { return S.Ground.CellSize; }, {0.05f, 10.f}},
        {[](FSceneSettings& S) -> float& { return S.Ground.MatchingHeight; }, {-2.f, 2.f}},
    };
    const FSceneSettings Before = Scene.GetSettings();
    const float Infinity = std::numeric_limits<float>::infinity();
    bool bRejected = true;
    bool bAccepted = true;
    for(const auto& Field : Fields)
    {
        const float Low = Field.second.first;
        const float High = Field.second.second;
        for(float Bad : {std::numeric_limits<float>::quiet_NaN(), Infinity, -Infinity, std::nextafter(Low, -Infinity), std::nextafter(High, Infinity)})
        {
            FSceneSettings Candidate = Before;
            // 同时携带一个合法改动，验证整体拒绝、不部分应用
            Candidate.CapsuleMaterials[3] = {0.7f, 0.3f};
            Candidate.Ground.bVisible = !Before.Ground.bVisible;
            Field.first(Candidate) = Bad;
            bRejected &= !Scene.SetSettings(Candidate) && SameSurfaceSettings(Scene.GetSettings(), Before);
        }
        for(float Good : {Low, High})
        {
            FSceneSettings Candidate = Before;
            Field.first(Candidate) = Good;
            bAccepted &= Scene.SetSettings(Candidate);
            FSceneSettings Current = Scene.GetSettings();
            bAccepted &= Field.first(Current) == Good;
            Scene.SetSettings(Before);
        }
    }
    Require(bRejected, "material/ground NaN, Inf and out-of-range rejected with whole config unchanged");
    Require(bAccepted, "material/ground boundary values r=.05/1 m=0/1 cell .05/10 height -2/2 accepted");
    Require(SameSurfaceSettings(Scene.GetSettings(), Before), "validation probe restores settings");
}

void ProbeMaterialsAndGround(FSandboxScene& Scene, FGpuRenderer& Renderer)
{
    FSceneCommand Restore;
    Restore.bRestoreDefaults = true;
    Restore.bReset = true;
    SetMode(Scene, ESceneMode::Physics);
    Scene.ApplyCommand(Restore);
    SetPaused(Scene, false);
    for(int Frame = 0; Frame <= 60; ++Frame) Scene.Tick({});
    SetPaused(Scene, true);
    const auto Frozen = Scene.GetSnapshot();
    const FSceneSettings Defaults = Scene.GetSettings();
    Require(Frozen.GroundPlane == std::array<double, 4>{0., 1., 0., 0.05} && Defaults.Ground.bVisible && Defaults.Ground.CellSize == 0.5f
            && Defaults.Ground.MatchingHeight == 0.f && SameMaterial(Defaults.Ground.Material, {0.f, 0.8f})
            && SameMaterial(Defaults.CharacterMaterial, {0.f, 0.1f}) && SameMaterial(Defaults.CapsuleMaterials[0], {0.f, 0.1f}),
            "defaults: physics plane = ContactPlane, ground visible 0.5m {0,0.8}, materials {0,0.1}");

    ProbeSettingsValidation(Scene);

    const FFramePair Default = Render(Scene, Renderer);
    Require(Parity(Default.Gpu, Default.Cpu), "default ground scene CPU/GPU parity");

    auto Configure = [&](bool bCharacter, bool bCapsules, bool bGround)
    {
        FSceneSettings Settings = Defaults;
        Settings.bShowCharacter = bCharacter;
        Settings.bShowCapsules = bCapsules;
        Settings.bShowContacts = false;
        Settings.Ground.bVisible = bGround;
        if(!Scene.SetSettings(Settings)) Require(false, "isolated object configuration accepted");
        return Settings;
    };
    const FSurfaceMaterial Changed = {1.f, 0.3f};

    // 胶囊：逐个改材质，变化应局限于单个胶囊且不改变其他对象设置
    Configure(false, false, false);
    const FFramePair Empty = Render(Scene, Renderer);
    FSceneSettings Capsules = Configure(false, true, false);
    const FFramePair CapsuleBase = Render(Scene, Renderer);
    const int CapsulePixels = Compare(Empty.Gpu, CapsuleBase.Gpu).Changed;
    int VisibleCapsules = 0;
    int BestCapsule = -1;
    int BestChanged = 0;
    bool bLocal = true;
    bool bOthersUntouched = true;
    for(int Capsule = 0; Capsule < SceneCapsuleCount; ++Capsule)
    {
        FSceneSettings Edited = Capsules;
        Edited.CapsuleMaterials[Capsule] = Changed;
        if(!Scene.SetSettings(Edited)) Require(false, "single capsule material accepted");
        for(int Other = 0; Other < SceneCapsuleCount; ++Other)
            if(Other != Capsule) bOthersUntouched &= SameMaterial(Scene.GetSettings().CapsuleMaterials[Other], Defaults.CapsuleMaterials[Other]);
        bOthersUntouched &= SameMaterial(Scene.GetSettings().CharacterMaterial, Defaults.CharacterMaterial)
            && SameMaterial(Scene.GetSettings().Ground.Material, Defaults.Ground.Material);
        FImageRGB8 Image;
        Scene.BuildRenderData(Renderer);
        Renderer.ReadbackRGB8(Image);
        const int ChangedPixels = Compare(CapsuleBase.Gpu, Image).Changed;
        VisibleCapsules += ChangedPixels > 0;
        bLocal &= ChangedPixels < CapsulePixels / 2;
        if(ChangedPixels > BestChanged) { BestChanged = ChangedPixels; BestCapsule = Capsule; }
    }
    std::cout << "capsule pixels=" << CapsulePixels << " capsules changed by own material=" << VisibleCapsules
              << " largest=" << BestCapsule << "(" << BestChanged << ")\n";
    Require(CapsulePixels > 0 && VisibleCapsules >= SceneCapsuleCount / 2 && bLocal && bOthersUntouched,
            "each capsule material affects only that capsule");
    FSceneSettings Edited = Capsules;
    Edited.CapsuleMaterials[BestCapsule] = Changed;
    Scene.SetSettings(Edited);
    const FFramePair CapsuleEdited = Render(Scene, Renderer);
    Require(Compare(CapsuleBase.Cpu, CapsuleEdited.Cpu).Changed > 0 && Parity(CapsuleEdited.Gpu, CapsuleEdited.Cpu),
            "capsule material edit CPU/GPU parity");
    Scene.SetSettings(Capsules);
    FFramePair Restored = Render(Scene, Renderer);
    Require(Restored.Gpu.Pixels == CapsuleBase.Gpu.Pixels && Restored.Cpu.Pixels == CapsuleBase.Cpu.Pixels, "capsule material restore reproduces image");

    // 角色与地面：各自独立材质
    struct FCase
    {
        const char* Name;
        bool bCharacter;
        bool bGround;
    };
    for(const FCase& Case : {FCase{"character", true, false}, FCase{"ground", false, true}})
    {
        FSceneSettings Base = Configure(Case.bCharacter, false, Case.bGround);
        const FFramePair Original = Render(Scene, Renderer);
        Require(Compare(Empty.Gpu, Original.Gpu).Changed > 0, Case.bGround ? "ground visible alone" : "character visible alone");
        FSceneSettings Material = Base;
        (Case.bGround ? Material.Ground.Material : Material.CharacterMaterial) = Changed;
        Scene.SetSettings(Material);
        const FFramePair Modified = Render(Scene, Renderer);
        bool bOthers = true;
        for(int Capsule = 0; Capsule < SceneCapsuleCount; ++Capsule)
            bOthers &= SameMaterial(Scene.GetSettings().CapsuleMaterials[Capsule], Defaults.CapsuleMaterials[Capsule]);
        bOthers &= SameMaterial(Case.bGround ? Scene.GetSettings().CharacterMaterial : Scene.GetSettings().Ground.Material,
                                Case.bGround ? Defaults.CharacterMaterial : Defaults.Ground.Material);
        Require(Compare(Original.Gpu, Modified.Gpu).Changed > 0 && Compare(Original.Cpu, Modified.Cpu).Changed > 0 && bOthers,
                Case.bGround ? "ground material changes ground only" : "character material changes character only");
        Require(Parity(Original.Gpu, Original.Cpu) && Parity(Modified.Gpu, Modified.Cpu),
                Case.bGround ? "ground material CPU/GPU parity" : "character material CPU/GPU parity");
        Scene.SetSettings(Base);
        Restored = Render(Scene, Renderer);
        Require(Restored.Gpu.Pixels == Original.Gpu.Pixels && Restored.Cpu.Pixels == Original.Cpu.Pixels,
                Case.bGround ? "ground material restore reproduces image" : "character material restore reproduces image");
    }

    // 格长改变画面且可恢复；关闭地面回到无地面画面
    FSceneSettings Cells = Configure(false, false, true);
    const FFramePair GroundOnly = Render(Scene, Renderer);
    Cells.Ground.CellSize = 1.5f;
    Scene.SetSettings(Cells);
    const FFramePair Coarse = Render(Scene, Renderer);
    Require(Compare(GroundOnly.Gpu, Coarse.Gpu).Changed > 0 && Parity(Coarse.Gpu, Coarse.Cpu), "cell size edit visible with CPU/GPU parity");
    Cells.Ground.bVisible = false;
    Scene.SetSettings(Cells);
    const FFramePair Hidden = Render(Scene, Renderer);
    Require(Hidden.Gpu.Pixels == Empty.Gpu.Pixels && Hidden.Cpu.Pixels == Empty.Cpu.Pixels, "hidden ground equals empty scene");

    const auto AfterEdits = Scene.GetSnapshot();
    Require(AfterEdits.NextFrame == Frozen.NextFrame && AfterEdits.RootPosition == Frozen.RootPosition && AfterEdits.bPaused,
            "paused material/ground edits do not advance skeleton");

    // Reset与模式切换保留设置；RestoreDefaults恢复全部
    FSceneSettings Custom = Defaults;
    Custom.CharacterMaterial = {0.4f, 0.6f};
    Custom.CapsuleMaterials[5] = {0.5f, 0.6f};
    Custom.Ground = {false, 0.3f, 1.f, {0.9f, 0.2f}};
    Require(Scene.SetSettings(Custom), "custom material/ground accepted");
    FSceneCommand Reset;
    Reset.bReset = true;
    Scene.ApplyCommand(Reset);
    Require(SameSurfaceSettings(Scene.GetSettings(), Custom), "reset keeps material/ground settings");
    SetMode(Scene, ESceneMode::MotionMatching);
    const auto Matching = Scene.GetSnapshot();
    Require(SameSurfaceSettings(Scene.GetSettings(), Custom) && Matching.GroundPlane == std::array<double, 4>{0., 1., 0., -static_cast<double>(0.3f)},
            "matching plane y=MatchingHeight with settings kept");
    SetPaused(Scene, false);
    for(int Step = 0; Step < 120; ++Step) Scene.Tick({0.f, -1.f, 0.f, 0.f});
    const auto Walked = Scene.GetSnapshot();
    Require(Walked.RootPosition != Matching.RootPosition && Walked.GroundPlane == Matching.GroundPlane, "root motion does not move ground");
    SetMode(Scene, ESceneMode::Physics);
    Require(SameSurfaceSettings(Scene.GetSettings(), Custom) && Scene.GetSnapshot().GroundPlane == Custom.ContactPlane,
            "physics plane restored after round trip with settings kept");
    Scene.ApplyCommand(Restore);
    Require(SameSurfaceSettings(Scene.GetSettings(), FSceneSettings{}), "restore defaults resets all materials and ground");

    // 只改地面显示/高度/格长/材质不影响运动
    SetPaused(Scene, false);
    auto Digests = [&]()
    {
        SetMode(Scene, ESceneMode::Physics);
        const uint64_t Physics = PhysicsDigest(Scene, 600);
        SetMode(Scene, ESceneMode::MotionMatching);
        const uint64_t Motion = MatchingDigest(Scene, 3000);
        return std::make_pair(Physics, Motion);
    };
    const auto DefaultDigest = Digests();
    FSceneSettings GroundEdited = Scene.GetSettings();
    GroundEdited.Ground = {false, 1.2f, 2.f, {1.f, 0.2f}};
    GroundEdited.CharacterMaterial = {1.f, 0.05f};
    Require(Scene.SetSettings(GroundEdited), "ground-only edit accepted");
    const auto EditedDigest = Digests();
    std::cout << "motion digest default physics600=" << std::hex << DefaultDigest.first << " matching3000=" << DefaultDigest.second
              << " ground-edited physics600=" << EditedDigest.first << " matching3000=" << EditedDigest.second << std::dec << '\n';
    Require(DefaultDigest == EditedDigest, "ground/material edits keep physics and matching motion digests");

    SetMode(Scene, ESceneMode::Physics);
    Scene.ApplyCommand(Restore);
    SetPaused(Scene, false);
}

void ProbeFrameCaptures(FSandboxScene& Scene, FGpuRenderer& Renderer, const std::string& Prefix)
{
    // 禁用地面截图与固定输入MM对拍由工具驱动，不添加产品CLI或改变输入路由
    FSceneCommand Restore;
    Restore.bRestoreDefaults = true;
    Restore.bReset = true;
    Scene.ApplyCommand(Restore);
    SetPaused(Scene, false);
    for(int Frame = 0; Frame <= 60; ++Frame) Scene.Tick({});
    auto Settings = Scene.GetSettings();
    Settings.Ground.bVisible = false;
    Require(Scene.SetSettings(Settings), "disable ground for corner capture");
    FFramePair Pair = Render(Scene, Renderer);
    Require(Parity(Pair.Gpu, Pair.Cpu), "ground disabled CPU/GPU parity");
    Save(Pair.Gpu, Prefix + "_hidden_gpu.png");
    Save(Pair.Cpu, Prefix + "_hidden_cpu.png");
    Scene.ApplyCommand(Restore);
    SetMode(Scene, ESceneMode::MotionMatching);
    for(int Frame = 0; Frame <= 300; ++Frame)
    {
        // 走120帧、停60帧、再横向转向；输入固定且沿用FSandboxScene::Tick
        const FSceneInput Input = Frame < 120 ? FSceneInput{0.f, -1.f, 0.f, 0.f}
            : Frame < 180 ? FSceneInput{} : FSceneInput{1.f, 0.f, 0.f, 0.f};
        Scene.Tick(Input);
        if(Frame != 60 && Frame != 300) continue;
        Pair = Render(Scene, Renderer);
        Require(Parity(Pair.Gpu, Pair.Cpu), "fixed-input MM ground CPU/GPU parity");
        const std::string Suffix = Frame == 60 ? "060.png" : "300.png";
        Save(Pair.Gpu, Prefix + "_mm_gpu_" + Suffix);
        Save(Pair.Cpu, Prefix + "_mm_cpu_" + Suffix);
    }
    SetMode(Scene, ESceneMode::Physics);
    Scene.ApplyCommand(Restore);
}

// ---- 逐脚触地标签文件 ----
std::vector<uint8_t> ReadBytes(const std::string& Path)
{
    std::ifstream File(Path, std::ios::binary);
    return std::vector<uint8_t>(std::istreambuf_iterator<char>(File), std::istreambuf_iterator<char>());
}
void WriteBytes(const std::string& Path, const std::vector<uint8_t>& Bytes)
{
    std::ofstream File(Path, std::ios::binary | std::ios::trunc);
    File.write(reinterpret_cast<const char*>(Bytes.data()), static_cast<std::streamsize>(Bytes.size()));
    if(!File) Require(false, "isolated corrupt copy written");
}
void Poke32(std::vector<uint8_t>& Bytes, size_t Offset, uint32_t Value)
{
    for(int Byte = 0; Byte < 4; ++Byte) Bytes[Offset + Byte] = static_cast<uint8_t>(Value >> (8 * Byte));
}

void ProbeFootContactFile(FSandboxScene& Scene, const std::string& Prefix)
{
    std::vector<uint8_t> Built;
    FFootContactStats Stats;
    std::string Report;
    Require(BuildFootContactsFile(Built, Stats, Report), "foot contacts rebuilt in memory");
    std::cout << "foot contacts " << Report << '\n';
    const std::vector<uint8_t> Stored = ReadBytes(GetSandboxResourcePath("foot_contacts.bin"));
    Require(Built == Stored, "in-memory labels byte-equal resources/foot_contacts.bin");
    std::vector<uint8_t> Flipped = Built;
    Flipped[FootContactHeaderSize + 1000] ^= FootContactLeft;
    Require(Flipped != Stored, "negative control single flipped label differs from stored file");
    bool bRates = true, bMirror = true;
    for(double Rate : Stats.ContactRate) bRates &= Rate >= 0.45 && Rate <= 0.60;
    for(double Agreement : Stats.MirrorAgreement) bMirror &= Agreement >= 0.999;
    Require(bRates, "per-foot contact rate in [0.45, 0.60]");
    Require(bMirror, "mirrored clip pairs left/right labels agree >= 99.9%");
    Require(Stats.StillFrames > 0 && Stats.StillBothRate >= 0.9, "root speed <0.1m/s frames double support >= 0.9");

    uint64_t SourceHash = 0;
    std::string Reason;
    Require(HashFileFnv1a64(GetSandboxResourcePath("long_motion.bin"), SourceHash, Reason), "long_motion.bin hashed");
    const int Frames = Stats.Frames;
    FFootContactFile Parsed;
    Require(ParseFootContacts(Built.data(), Built.size(), Frames, SourceHash, Parsed, Reason) == EFootContactError::None
            && static_cast<int>(Parsed.Labels.size()) == Frames, "valid label file parses");
    Require(ParseFootContacts(Built.data(), Built.size(), Frames - 1, SourceHash, Parsed, Reason) == EFootContactError::FrameCount
            && ParseFootContacts(Built.data(), Built.size(), Frames, SourceHash ^ 1u, Parsed, Reason) == EFootContactError::Hash,
            "caller frame count and source hash mismatches rejected");

    // 损坏副本：内存解析与隔离临时文件加载均须给出对应错误类别；临时文件用后删除。
    struct FCorrupt
    {
        const char* Name;
        EFootContactError Expected;
        std::vector<uint8_t> Bytes;
    };
    std::vector<FCorrupt> Cases;
    auto Variant = [&](const char* Name, EFootContactError Expected, const std::function<void(std::vector<uint8_t>&)>& Edit)
    {
        std::vector<uint8_t> Bytes = Built;
        Edit(Bytes);
        Cases.push_back({Name, Expected, std::move(Bytes)});
    };
    Variant("magic", EFootContactError::Magic, [](std::vector<uint8_t>& Bytes) { Bytes[0] ^= 1u; });
    Variant("version", EFootContactError::Version, [](std::vector<uint8_t>& Bytes) { Poke32(Bytes, 8, FootContactVersion + 1); });
    Variant("truncated", EFootContactError::Truncated, [](std::vector<uint8_t>& Bytes) { Bytes.pop_back(); });
    Variant("header", EFootContactError::Truncated, [](std::vector<uint8_t>& Bytes) { Bytes.resize(40); });
    Variant("frames", EFootContactError::FrameCount, [Frames](std::vector<uint8_t>& Bytes)
            { Poke32(Bytes, 12, static_cast<uint32_t>(Frames - 1)); Bytes.pop_back(); });
    Variant("hash", EFootContactError::Hash, [](std::vector<uint8_t>& Bytes) { Bytes[56] ^= 1u; });
    Variant("label", EFootContactError::Label, [](std::vector<uint8_t>& Bytes) { Bytes[FootContactHeaderSize + 7] = 4; });
    Variant("trailing", EFootContactError::Trailing, [](std::vector<uint8_t>& Bytes) { Bytes.push_back(0); });
    Variant("bones", EFootContactError::Bones, [](std::vector<uint8_t>& Bytes) { Poke32(Bytes, 16, 3); });
    bool bMemory = true, bFiles = true;
    for(const FCorrupt& Case : Cases)
    {
        const EFootContactError Memory = ParseFootContacts(Case.Bytes.data(), Case.Bytes.size(), Frames, SourceHash, Parsed, Reason);
        const std::string Path = Prefix + "_corrupt_" + Case.Name + ".bin";
        WriteBytes(Path, Case.Bytes);
        const EFootContactError Loaded = LoadFootContacts(Path, Frames, SourceHash, Parsed, Reason);
        std::filesystem::remove(Path);
        std::cout << "  corrupt " << Case.Name << " memory=" << GetFootContactErrorName(Memory) << " file="
                  << GetFootContactErrorName(Loaded) << " (" << Reason << ")\n";
        bMemory &= Memory == Case.Expected;
        bFiles &= Loaded == Case.Expected && !std::filesystem::exists(Path);
    }
    Require(bMemory, "corrupt labels (magic/version/truncated/header/frames/hash/label/trailing/bones) rejected in memory");
    Require(bFiles, "isolated corrupt label files rejected by loader and removed");
    const std::string Missing = Prefix + "_missing_foot_contacts.bin";
    Require(!std::filesystem::exists(Missing) && LoadFootContacts(Missing, Frames, SourceHash, Parsed, Reason) == EFootContactError::Open,
            "missing label file rejected as open failure");

    // 标签缺失或损坏：MM不可用并给出原因，Physics照常且摘要与正常场景一致。
    FSceneCommand Restore;
    Restore.bRestoreDefaults = true;
    Restore.bReset = true;
    SetMode(Scene, ESceneMode::Physics);
    Scene.ApplyCommand(Restore);
    SetPaused(Scene, false);
    const uint64_t Physics = PhysicsDigest(Scene, 600);
    const std::string Corrupt = Prefix + "_corrupt_scene.bin";
    const auto HashCase = std::find_if(Cases.begin(), Cases.end(), [](const FCorrupt& Case) { return Case.Expected == EFootContactError::Hash; });
    WriteBytes(Corrupt, HashCase->Bytes);
    bool bUnavailable = true;
    for(const std::string& Path : {Missing, Corrupt})
    {
        FSandboxScene Broken(Path);
        const bool bAvailable = Broken.GetSnapshot().bMatchingAvailable;
        FSceneCommand Enter;
        Enter.bSetMode = true;
        Enter.Mode = ESceneMode::MotionMatching;
        const bool bEntered = Broken.ApplyCommand(Enter);
        const auto After = Broken.GetSnapshot();
        const uint64_t BrokenPhysics = PhysicsDigest(Broken, 600);
        std::cout << "  broken labels message: " << After.Message << " physics600=" << std::hex << BrokenPhysics << std::dec << '\n';
        bUnavailable &= !bAvailable && !bEntered && After.Mode == ESceneMode::Physics
            && After.Message.find("逐脚触地标签无效") != std::string::npos && BrokenPhysics == Physics;
    }
    std::filesystem::remove(Corrupt);
    std::cout << "foot label physics600=" << std::hex << Physics << std::dec << '\n';
    Require(bUnavailable && !std::filesystem::exists(Corrupt),
            "missing/corrupt labels: matching unavailable with reason, physics digest unchanged");
    SetMode(Scene, ESceneMode::MotionMatching);
}

// ---- 稳定两骨IK受控链（独立锚点） ----
using FVec3d = std::array<double, 3>;
using FQuatd = std::array<double, 4>;
FVec3d SubD(const FVec3d& A, const FVec3d& B) { return {A[0] - B[0], A[1] - B[1], A[2] - B[2]}; }
FVec3d AddD(const FVec3d& A, const FVec3d& B) { return {A[0] + B[0], A[1] + B[1], A[2] + B[2]}; }
FVec3d ScaleD(const FVec3d& A, double S) { return {A[0] * S, A[1] * S, A[2] * S}; }
FVec3d CrossD(const FVec3d& A, const FVec3d& B)
{
    return {A[1] * B[2] - A[2] * B[1], A[2] * B[0] - A[0] * B[2], A[0] * B[1] - A[1] * B[0]};
}
double NormD(const FVec3d& A) { return std::sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2]); }
double DistanceD(const FVec3d& A, const FVec3d& B) { return NormD(SubD(A, B)); }
// q(w,x,y,z)作用于v：v + 2w(u x v) + 2u x (u x v)
FVec3d RotateD(const FQuatd& Q, const FVec3d& V)
{
    const FVec3d U = {Q[1], Q[2], Q[3]};
    const FVec3d T = ScaleD(CrossD(U, V), 2.);
    return AddD(AddD(V, ScaleD(T, Q[0])), CrossD(U, T));
}
bool UnitD(const FQuatd& Q)
{
    const double Norm = std::sqrt(Q[0] * Q[0] + Q[1] * Q[1] + Q[2] * Q[2] + Q[3] * Q[3]);
    return std::isfinite(Norm) && std::fabs(Norm - 1.) <= 1e-9;
}
// 输出旋转单位、按输入偏移前向运动学复现Mid/End、骨长保持（<=1e-5m）。
bool ConsistentIk(const FTwoBoneIkInput& In, const FTwoBoneIkOutput& Out)
{
    const FVec3d Mid = AddD(In.Root, RotateD(Out.UpperRotation, In.MidOffset));
    const FVec3d End = AddD(Out.Mid, RotateD(Out.LowerRotation, In.EndOffset));
    return UnitD(Out.UpperRotation) && UnitD(Out.LowerRotation) && DistanceD(Mid, Out.Mid) <= 1e-9 && DistanceD(End, Out.End) <= 1e-9
        && std::fabs(DistanceD(Out.Mid, In.Root) - NormD(In.MidOffset)) <= 1e-5
        && std::fabs(DistanceD(Out.End, Out.Mid) - NormD(In.EndOffset)) <= 1e-5 && std::isfinite(Out.Unreachable);
}

void ProbeTwoBoneIk()
{
    // 锚点1：两段0.5m直腿竖直向下，目标(0,-0.6,0)，参考轴+x：内角cos=0.6，膝偏0.4、沿线0.3，
    // 弯曲法向(c-a)x(b-a)=(0.24,0,0)与参考轴同向，故膝在(0,-0.3,-0.4)。
    FTwoBoneIkInput Straight;
    Straight.MidOffset = {0., -0.5, 0.};
    Straight.EndOffset = {0., -0.5, 0.};
    Straight.Target = {0., -0.6, 0.};
    Straight.ReferenceAxis = {1., 0., 0.};
    FTwoBoneIkOutput Bent;
    Require(SolveTwoBoneIk(Straight, Bent) && Bent.bDegenerateAxis && DistanceD(Bent.Mid, {0., -0.3, -0.4}) <= 1e-5
            && DistanceD(Bent.End, Straight.Target) <= 1e-5 && Bent.Unreachable <= 1e-12 && ConsistentIk(Straight, Bent),
            "IK anchor: straight leg bends about reference axis, knee (0,-0.3,-0.4)");
    // 锚点2：以锚点1输出为输入、参考轴取反，膝翻到另一侧(0,-0.3,+0.4)；参考轴不变时原样保持。
    FTwoBoneIkInput Flip = Straight;
    Flip.UpperRotation = Bent.UpperRotation;
    Flip.LowerRotation = Bent.LowerRotation;
    FTwoBoneIkOutput Kept;
    Require(SolveTwoBoneIk(Flip, Kept) && !Kept.bDegenerateAxis && DistanceD(Kept.Mid, Bent.Mid) <= 1e-5
            && DistanceD(Kept.End, Flip.Target) <= 1e-5 && ConsistentIk(Flip, Kept), "IK anchor: satisfied target is a fixed point");
    Flip.ReferenceAxis = {-1., 0., 0.};
    FTwoBoneIkOutput Flipped;
    Require(SolveTwoBoneIk(Flip, Flipped) && !Flipped.bDegenerateAxis && DistanceD(Flipped.Mid, {0., -0.3, 0.4}) <= 1e-5
            && DistanceD(Flipped.End, Flip.Target) <= 1e-5 && ConsistentIk(Flip, Flipped),
            "IK anchor: reversed reference axis flips knee to (0,-0.3,+0.4)");
    // 不可达：直链当前伸展1.0，D=1.0；弯链当前伸展0.6，D=0.999。
    FTwoBoneIkInput Far = Straight;
    Far.Target = {0., -2., 0.};
    FTwoBoneIkOutput Clamped;
    Require(SolveTwoBoneIk(Far, Clamped) && DistanceD(Clamped.End, {0., -1., 0.}) <= 1e-5 && std::fabs(Clamped.Unreachable - 1.) <= 1e-9
            && ConsistentIk(Far, Clamped), "IK anchor: unreachable straight leg clamped to 1.0, unreachable 1.0");
    Far.UpperRotation = Bent.UpperRotation;
    Far.LowerRotation = Bent.LowerRotation;
    Require(SolveTwoBoneIk(Far, Clamped) && DistanceD(Clamped.End, {0., -0.999, 0.}) <= 1e-5
            && std::fabs(Clamped.Unreachable - 1.001) <= 1e-9 && ConsistentIk(Far, Clamped),
            "IK anchor: unreachable bent leg clamped to 0.999(l1+l2), unreachable 1.001");
    // 退化：目标在根点（钳到最小距离1e-3并保持方向）、目标在正上方（180度分支）、非有限与零骨长拒绝。
    FTwoBoneIkInput AtRoot = Straight;
    AtRoot.Target = {0., 0., 0.};
    FTwoBoneIkOutput Degenerate;
    Require(SolveTwoBoneIk(AtRoot, Degenerate) && std::fabs(DistanceD(Degenerate.End, AtRoot.Root) - 1e-3) <= 1e-7
            && ConsistentIk(AtRoot, Degenerate), "IK degenerate: target at root finite, unit, reach clamped to 1e-3");
    FTwoBoneIkInput Above = Straight;
    Above.Target = {0., 0.6, 0.};
    Require(SolveTwoBoneIk(Above, Degenerate) && DistanceD(Degenerate.End, Above.Target) <= 1e-5 && ConsistentIk(Above, Degenerate),
            "IK degenerate: target opposite current reach (180deg swing) reached");
    FTwoBoneIkInput Invalid = Straight;
    Invalid.Target[0] = std::numeric_limits<double>::quiet_NaN();
    FTwoBoneIkInput ZeroBone = Straight;
    ZeroBone.MidOffset = {0., 0., 0.};
    Require(!SolveTwoBoneIk(Invalid, Degenerate) && !SolveTwoBoneIk(ZeroBone, Degenerate), "IK rejects nonfinite target and zero bone");

    // 随机受控链：可达目标末端误差<=1e-4m；不可达目标沿目标方向钳到D且Unreachable=|t-a|-D；直腿随机朝向退化。
    uint32_t Seed = 917u;
    auto Uniform = [&](double Low, double High)
    {
        Seed = Seed * 1664525u + 1013904223u;
        return Low + (High - Low) * ((Seed >> 8) / 16777216.);
    };
    auto RandomUnit = [&]()
    {
        FVec3d V;
        do V = {Uniform(-1., 1.), Uniform(-1., 1.), Uniform(-1., 1.)}; while(NormD(V) < 0.1 || NormD(V) > 1.);
        return ScaleD(V, 1. / NormD(V));
    };
    auto RandomRotation = [&]()
    {
        FQuatd Q;
        double Norm = 0.;
        do { Q = {Uniform(-1., 1.), Uniform(-1., 1.), Uniform(-1., 1.), Uniform(-1., 1.)}; Norm = std::sqrt(Q[0] * Q[0] + Q[1] * Q[1] + Q[2] * Q[2] + Q[3] * Q[3]); }
        while(Norm < 0.1 || Norm > 1.);
        for(double& Value : Q) Value /= Norm;
        return Q;
    };
    double MaxReachError = 0., MaxClampError = 0., MaxStraightError = 0.;
    bool bReachable = true, bUnreachable = true, bStraight = true;
    for(int Case = 0; Case < 3000; ++Case)
    {
        const int Kind = Case % 3;
        FTwoBoneIkInput In;
        In.Root = {Uniform(-5., 5.), Uniform(0., 2.), Uniform(-5., 5.)};
        In.UpperRotation = RandomRotation();
        In.LowerRotation = Kind == 2 ? In.UpperRotation : RandomRotation();
        const FVec3d Direction = RandomUnit();
        In.MidOffset = ScaleD(Direction, Uniform(0.2, 0.6));
        In.EndOffset = ScaleD(Kind == 2 ? Direction : RandomUnit(), Uniform(0.2, 0.6));
        In.ReferenceAxis = RandomUnit();
        const double L1 = NormD(In.MidOffset), L2 = NormD(In.EndOffset);
        const double MinReach = std::fabs(L1 - L2) + 1e-3 * (L1 + L2);
        const double Raw = Kind == 1 ? Uniform(1.01 * (L1 + L2), 3. * (L1 + L2)) : Uniform(MinReach + 0.01 * (L1 + L2), 0.99 * (L1 + L2));
        In.Target = AddD(In.Root, ScaleD(RandomUnit(), Raw));
        FTwoBoneIkOutput Out;
        const bool bSolved = SolveTwoBoneIk(In, Out) && ConsistentIk(In, Out);
        if(Kind == 1)
        {
            const FVec3d Mid = AddD(In.Root, RotateD(In.UpperRotation, In.MidOffset));
            const double Current = DistanceD(AddD(Mid, RotateD(In.LowerRotation, In.EndOffset)), In.Root);
            const double D = std::min(L1 + L2, std::max(0.999 * (L1 + L2), Current));
            const FVec3d Expected = AddD(In.Root, ScaleD(SubD(In.Target, In.Root), D / Raw));
            const double Error = DistanceD(Out.End, Expected);
            MaxClampError = std::max(MaxClampError, Error);
            bUnreachable &= bSolved && Error <= 1e-4 && std::fabs(Out.Unreachable - (Raw - D)) <= 1e-9;
        }
        else
        {
            const double Error = DistanceD(Out.End, In.Target);
            double& MaxError = Kind == 0 ? MaxReachError : MaxStraightError;
            bool& bPass = Kind == 0 ? bReachable : bStraight;
            MaxError = std::max(MaxError, Error);
            bPass &= bSolved && Error <= 1e-4 && Out.Unreachable <= 1e-12 && (Kind == 0 || Out.bDegenerateAxis);
        }
    }
    std::cout << "two-bone IK random maxReachError=" << MaxReachError << " maxStraightError=" << MaxStraightError
              << " maxClampError=" << MaxClampError << '\n';
    Require(bReachable, "IK random reachable targets: end error <=1e-4m, bone length <=1e-5m, unit rotations");
    Require(bStraight, "IK random straight legs: degenerate axis fallback reaches target with unit rotations");
    Require(bUnreachable, "IK random unreachable targets: clamped along target direction, no NaN");
}

// ---- 脚尖支点IK：闭式锚点、soft连续性与边界负对照 ----
bool ConsistentPivot(const FFootPivotIkInput& In, const FFootPivotIkOutput& Out)
{
    return ConsistentIk(In.Chain, Out.Chain) && UnitD(Out.FootRotation)
        && DistanceD(AddD(Out.Chain.End, RotateD(Out.FootRotation, In.ToeOffset)), Out.Toe) <= 1e-9
        && std::fabs(DistanceD(Out.Chain.End, Out.Toe) - NormD(In.ToeOffset)) <= 1e-5
        && Out.LiftAngle >= 0. && Out.LiftAngle <= 3.14159265358979323846 / 4. + 1e-12 && std::isfinite(Out.Residual)
        && std::fabs(Out.Residual - DistanceD(Out.Toe, In.ToeTarget)) <= 1e-9
        && Out.Chain.End[1] - Out.Toe[1] >= -RotateD(In.FootRotation, In.ToeOffset)[1] - 1e-9;
}
void ProbeFootPivotIk()
{
    constexpr double PiD = 3.14159265358979323846;
    FFootPivotIkInput In;
    In.Chain.Root = {0., 0.8, 0.45};
    In.Chain.MidOffset = {0., -0.4, -std::sqrt(0.45 * 0.45 - 0.4 * 0.4)};
    In.Chain.EndOffset = {0., -0.4, std::sqrt(0.45 * 0.45 - 0.4 * 0.4)};
    In.ToeOffset = {0., -0.075, std::sqrt(0.15 * 0.15 - 0.075 * 0.075)};
    const auto AnimatedToe = AddD(AddD(AddD(In.Chain.Root, In.Chain.MidOffset), In.Chain.EndOffset), In.ToeOffset);
    In.ToeTarget = AnimatedToe;
    FFootPivotIkOutput Out;
    Require(SolveFootPivotIk(In, Out) && ConsistentPivot(In, Out) && Out.LiftAngle <= 1e-12
            && DistanceD(Out.Toe, AnimatedToe) <= 1e-9 && !Out.bConflict, "pivot anchor: animation target is a no-lift fixed point");

    // 两段0.45m、脚掌0.15m、髋(0,.8,.45)、锁点0；膝角带上端取D=.88。
    // 用独立二维余弦闭式算抬跟角：ψ - acos((|H|²+|u|²-D²)/(2|H||u|))。
    In.ToeTarget = {0., 0., 0.};
    In.PreviousKneeAngle = std::acos((2. * 0.45 * 0.45 - 0.88 * 0.88) / (2. * 0.45 * 0.45)) - 14.5 * PiD / 180.;
    const double Psi = std::atan2(0.075 * 0.45 + In.ToeOffset[2] * 0.8, 0.075 * 0.8 - In.ToeOffset[2] * 0.45);
    const double LiftExpected = Psi - std::acos((0.8 * 0.8 + 0.45 * 0.45 + 0.15 * 0.15 - 0.88 * 0.88)
                                               / (2. * 0.15 * std::hypot(0.8, 0.45)));
    Require(SolveFootPivotIk(In, Out) && ConsistentPivot(In, Out) && !Out.bConflict && Out.Residual <= 1e-9
            && std::fabs(Out.LiftAngle - LiftExpected) <= 1e-9 && std::fabs(DistanceD(Out.Chain.End, In.Chain.Root) - 0.88) <= 1e-9,
            "pivot anchor: toe fixed, reach .88m, lift matches independent closed form (~18.57deg)");
    auto Bad = Out;
    Bad.Toe[0] += 0.002;
    Require(!ConsistentPivot(In, Bad), "pivot negative control: rewriting toe without FK rejected");
    Bad = Out;
    Bad.FootRotation = {1., 0., 0., 0.};
    Require(!ConsistentPivot(In, Bad), "pivot negative control: mismatched foot rotation rejected");

    auto Landing = In;
    Landing.Chain.Root = {0., 0.8, -0.7};
    Landing.PreviousKneeAngle = PiD;
    Require(SolveFootPivotIk(Landing, Out) && ConsistentPivot(Landing, Out) && Out.bConflict && Out.LiftAngle <= 1e-12
            && std::fabs(Out.Residual - (std::hypot(0.8 - 0.075, -0.7 + Landing.ToeOffset[2]) - 0.8991)) <= 1e-9,
            "pivot boundary: 20deg conditioning margin limits landing lift; knee-priority residual reported");
    const double Epsilon = 1e-6, Start = 0.9 - 0.05;
    const double Left = (SoftenFootReach(Start, 0.9) - SoftenFootReach(Start - Epsilon, 0.9)) / Epsilon;
    const double Right = (SoftenFootReach(Start + Epsilon, 0.9) - SoftenFootReach(Start, 0.9)) / Epsilon;
    Require(std::fabs(SoftenFootReach(Start, 0.9) - Start) <= 1e-12 && std::fabs(Left - Right) <= 2e-5,
            "soft anchor: value and first derivative continuous at 5cm zone boundary");
    bool bSoft = true, bPivot = true;
    int Conflicts = 0, Exact = 0;
    double PreviousSoft = SoftenFootReach(0., 0.9), PreviousKnee = -1.;
    for(int Index = 1; Index <= 3000; ++Index)
    {
        const double D = Index * 0.0005;
        const double Soft = SoftenFootReach(D, 0.9);
        bSoft &= Soft >= PreviousSoft && Soft - PreviousSoft <= 0.0005 + 1e-12 && Soft <= 0.9;
        PreviousSoft = Soft;
        auto Sample = In;
        Sample.Chain.Root = {0.04 * std::sin(Index * 0.03), 0.8, 0.8 * std::sin(Index * 0.02)};
        Sample.PreviousKneeAngle = PreviousKnee;
        Sample.SoftWeight = (Index % 101) / 100.; // 覆盖淡出的完整soft权重区间。
        const bool bSolved = SolveFootPivotIk(Sample, Out) && ConsistentPivot(Sample, Out);
        const auto A = SubD(Sample.Chain.Root, Out.Chain.Mid), B = SubD(Out.Chain.End, Out.Chain.Mid);
        const double MeasuredKnee = std::atan2(NormD(CrossD(A, B)), A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
        bPivot &= bSolved && (Out.bConflict || Out.Residual <= 1e-9) && std::fabs(MeasuredKnee - Out.KneeAngle) <= 1e-9
            && (PreviousKnee < 0. || std::fabs(MeasuredKnee - PreviousKnee) <= 14.5 * PiD / 180. + 1e-9);
        if(bSolved) { Conflicts += Out.bConflict; Exact += !Out.bConflict; PreviousKnee = MeasuredKnee; }
    }
    Require(bSoft, "soft scan: monotonic, 1-Lipschitz, bounded reach (3000 samples)");
    Require(bPivot && Conflicts > 0 && Exact > 0, "pivot moving targets: knee step <=14.5deg, exact/conflict paths both covered, FK/unit/length valid");
    auto Vertical = In;
    Vertical.ToeOffset = {0., -0.15, 0.};
    Vertical.PreviousLiftAxis = {0., 0., 0.};
    Require(SolveFootPivotIk(Vertical, Out) && ConsistentPivot(Vertical, Out), "pivot boundary: vertical foot with degenerate previous axis finite");
    auto AtHip = In;
    AtHip.ToeTarget = AddD(AtHip.Chain.Root, AtHip.ToeOffset);
    Require(SolveFootPivotIk(AtHip, Out) && ConsistentPivot(AtHip, Out), "pivot boundary: requested ankle at hip finite");
    auto Invalid = In;
    Invalid.ToeTarget[0] = std::numeric_limits<double>::quiet_NaN();
    auto Zero = In;
    Zero.ToeOffset = {};
    auto ZeroLeg = In;
    ZeroLeg.Chain.MidOffset = {};
    Require(!SolveFootPivotIk(Invalid, Out) && !SolveFootPivotIk(Zero, Out) && !SolveFootPivotIk(ZeroLeg, Out)
            && std::isnan(SoftenFootReach(-1., .9)), "pivot boundary: NaN, zero foot/leg and invalid soft distance rejected");
}

// ---- 脚锁IK场景对照（IK关/开、同一输入脚本） ----
struct FFootStep
{
    std::array<float, 3> Root = {};
    float Yaw = 0.f;
    int Source = -1;
    int SearchCount = 0;
    int Query = -1;
    int Result = -1;
    int Contacts = 0;
    bool bFallback = false;
    bool bBlending = false;
    std::array<FFootView, 2> Feet;
    std::array<FBoneView, SceneCapsuleCount> Bones;
    double IkMs = 0.;
};
FFootStep RecordFoot(const FSceneSnapshot& Snapshot)
{
    FFootStep Step;
    Step.Root = Snapshot.RootPosition;
    Step.Yaw = Snapshot.RootYaw;
    Step.Source = Snapshot.SourceFrame;
    Step.SearchCount = Snapshot.SearchCount;
    Step.Query = Snapshot.SearchQueryContacts;
    Step.Result = Snapshot.SearchResultContacts;
    Step.Contacts = Snapshot.FootContacts;
    Step.bFallback = Snapshot.bSearchContactFallback;
    Step.bBlending = Snapshot.bBlending;
    Step.Feet = Snapshot.Feet;
    Step.Bones = Snapshot.Bones;
    Step.IkMs = Snapshot.FootIkMs;
    return Step;
}
// 逐字段摘要，不哈希结构体填充或计时；用于保存本次调整前的关闭路径和检索序列锚点。
uint64_t FootRunDigest(const std::vector<FFootStep>& Run, bool bPose)
{
    FMotionDigest Digest;
    for(const FFootStep& Step : Run)
    {
        Digest.Bytes(Step.Root.data(), sizeof(Step.Root));
        Digest.Bytes(&Step.Yaw, sizeof(Step.Yaw));
        for(int Value : {Step.Source, Step.SearchCount, Step.Query, Step.Result, Step.Contacts}) Digest.Bytes(&Value, sizeof(Value));
        if(!bPose) continue;
        for(const auto& Bone : Step.Bones)
        {
            Digest.Bytes(Bone.Position.data(), sizeof(Bone.Position));
            Digest.Bytes(Bone.Rotation.data(), sizeof(Bone.Rotation));
        }
        for(const auto& Foot : Step.Feet)
        {
            for(const auto& Position : {Foot.HipPosition, Foot.KneePosition, Foot.AnklePosition, Foot.ToePosition})
                Digest.Bytes(Position.data(), sizeof(Position));
        }
    }
    return Digest.Hash;
}
using FInputScript = std::function<FSceneInput(int)>;
// 固定脚本（相机朝向不变，输入即相机局部方向）：起步直行、90度右转、180度折返、停步、再起步后退、停步。
FInputScript FixedFootScript()
{
    return [](int Step) -> FSceneInput
    {
        if(Step < 120) return {0.f, -1.f, 0.f, 0.f};
        if(Step < 210) return {1.f, 0.f, 0.f, 0.f};
        if(Step < 330) return {-1.f, 0.f, 0.f, 0.f};
        if(Step < 450) return {};
        if(Step < 540) return {0.f, 1.f, 0.f, 0.f};
        return {};
    };
}
// 与RunRandom同一生成器：分段常量八方向/停止/模拟量/偶发非有限值与相机环绕。
FInputScript RandomFootScript(uint32_t Seed)
{
    struct FState
    {
        uint32_t Seed = 0;
        int Remaining = 0;
        FSceneInput Input;
    };
    auto State = std::make_shared<FState>();
    State->Seed = Seed;
    return [State](int) -> FSceneInput
    {
        auto Next = [&]() { State->Seed = State->Seed * 1664525u + 1013904223u; return State->Seed >> 8; };
        if(State->Remaining-- <= 0)
        {
            State->Remaining = 10 + Next() % 150;
            const uint32_t Kind = Next() % 16;
            FSceneInput& Input = State->Input;
            Input = {};
            if(Kind < 9)
            {
                static const float Dirs[9][2] = {{0, -1}, {1, -1}, {1, 0}, {1, 1}, {0, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {0, 0}};
                Input.MoveX = Dirs[Kind][0];
                Input.MoveY = Dirs[Kind][1];
            }
            else if(Kind < 14)
            {
                Input.MoveX = (Next() % 2001) / 1000.f - 1.f;
                Input.MoveY = (Next() % 2001) / 1000.f - 1.f;
            }
            else if(Kind == 14) Input.MoveX = Input.MoveY = std::numeric_limits<float>::quiet_NaN();
            if(Next() % 4 == 0) Input.LookX = (Next() % 2001) / 1000.f - 1.f;
        }
        return State->Input;
    };
}
std::vector<FFootStep> RunFootScript(FSandboxScene& Scene, bool bFootIk, int Steps, const FInputScript& Script)
{
    FSceneSettings Settings = Scene.GetSettings();
    Settings.bFootIk = bFootIk;
    if(!Scene.SetSettings(Settings)) Require(false, "foot IK setting accepted");
    Scene.Reset();
    std::vector<FFootStep> Run;
    Run.reserve(static_cast<size_t>(Steps) + 1);
    Run.push_back(RecordFoot(Scene.GetSnapshot()));
    for(int Step = 0; Step < Steps; ++Step)
    {
        Scene.Tick(Script(Step));
        Run.push_back(RecordFoot(Scene.GetSnapshot()));
    }
    return Run;
}
float HorizontalDistance(const std::array<float, 3>& A, const std::array<float, 3>& B) { return std::hypot(A[0] - B[0], A[2] - B[2]); }
float Distance3(const std::array<float, 3>& A, const std::array<float, 3>& B)
{
    const float X = A[0] - B[0], Y = A[1] - B[1], Z = A[2] - B[2];
    return std::sqrt(X * X + Y * Y + Z * Z);
}
// 膝内角（髋-膝-踝），度。
float KneeDegrees(const FFootView& Foot)
{
    const FVec3d U = SubD({Foot.HipPosition[0], Foot.HipPosition[1], Foot.HipPosition[2]},
                          {Foot.KneePosition[0], Foot.KneePosition[1], Foot.KneePosition[2]});
    const FVec3d V = SubD({Foot.AnklePosition[0], Foot.AnklePosition[1], Foot.AnklePosition[2]},
                          {Foot.KneePosition[0], Foot.KneePosition[1], Foot.KneePosition[2]});
    return Degrees(static_cast<float>(std::atan2(NormD(CrossD(U, V)), U[0] * V[0] + U[1] * V[1] + U[2] * V[2])));
}
struct FFootRunMetrics
{
    // 标签前后两步均触地时显示脚尖的累计水平位移（两脚合计）
    double LabelSlip = 0.;
    float MaxKneeDelta = 0.f;
    // 锁定持续步（前后两步均锁定且锁点不变）
    int LockedSteps = 0;
    int LockedOver = 0;
    float MaxLockedStep = 0.f;
    float MaxLockOffset = 0.f;
    float MaxLockHeightError = 0.f;
    // 显示脚尖3D步长 - 动画脚尖3D步长
    int ExtraOver = 0;
    float MaxExtra = 0.f;
    int FadeSteps = 0;
    float MaxFadeExtra = 0.f;
    // 分类诊断：锁定超限步中含几何冲突(bConflict)的步数；膝角步变按 锁定/淡出/状态切换/自由 归类
    int LockedOverClamped = 0;
    float MaxLockedStepUnclamped = 0.f;
    std::array<float, 4> MaxKneeByState = {};
    std::array<float, 4> MaxClampByKnee = {};
    // 标签触地滑移中前后两步均锁定的部分
    double LabelSlipLocked = 0.;
};
FFootRunMetrics MeasureFootRun(const std::vector<FFootStep>& Run)
{
    FFootRunMetrics Metrics;
    for(size_t Index = 1; Index < Run.size(); ++Index)
    {
        for(int Side = 0; Side < 2; ++Side)
        {
            const FFootView& Previous = Run[Index - 1].Feet[Side];
            const FFootView& Current = Run[Index].Feet[Side];
            const float Horizontal = HorizontalDistance(Current.ToePosition, Previous.ToePosition);
            if(Previous.bContactLabel && Current.bContactLabel)
            {
                Metrics.LabelSlip += Horizontal;
                if(Previous.bLocked && Current.bLocked) Metrics.LabelSlipLocked += Horizontal;
            }
            const float KneeDelta = std::fabs(KneeDegrees(Current) - KneeDegrees(Previous));
            Metrics.MaxKneeDelta = std::max(Metrics.MaxKneeDelta, KneeDelta);
            const int State = Previous.bLocked != Current.bLocked ? 2 : Current.bLocked ? 0 : (Previous.bFading || Current.bFading) ? 1 : 3;
            if(KneeDelta > Metrics.MaxKneeByState[State])
            {
                Metrics.MaxKneeByState[State] = KneeDelta;
                Metrics.MaxClampByKnee[State] = std::max(Previous.ClampError, Current.ClampError);
            }
            if(Current.bLocked)
            {
                Metrics.MaxLockOffset = std::max(Metrics.MaxLockOffset, HorizontalDistance(Current.ToePosition, Current.LockPosition));
                Metrics.MaxLockHeightError = std::max(Metrics.MaxLockHeightError, std::fabs(Current.ToePosition[1] - Current.AnimToePosition[1]));
            }
            if(Previous.bLocked && Current.bLocked && Previous.LockPosition == Current.LockPosition)
            {
                ++Metrics.LockedSteps;
                Metrics.LockedOver += Horizontal > 1e-3f;
                Metrics.MaxLockedStep = std::max(Metrics.MaxLockedStep, Horizontal);
                const bool bClamped = Previous.bConflict || Current.bConflict;
                Metrics.LockedOverClamped += Horizontal > 1e-3f && bClamped;
                if(!bClamped) Metrics.MaxLockedStepUnclamped = std::max(Metrics.MaxLockedStepUnclamped, Horizontal);
            }
            const float Extra = Distance3(Current.ToePosition, Previous.ToePosition)
                - Distance3(Current.AnimToePosition, Previous.AnimToePosition);
            Metrics.MaxExtra = std::max(Metrics.MaxExtra, Extra);
            Metrics.ExtraOver += Extra > 0.02f;
            if(Previous.bFading || Current.bFading)
            {
                ++Metrics.FadeSteps;
                Metrics.MaxFadeExtra = std::max(Metrics.MaxFadeExtra, Extra);
            }
        }
    }
    return Metrics;
}

bool SameBone(const FBoneView& A, const FBoneView& B) { return A.Position == B.Position && A.Rotation == B.Rotation; }
FVec3d ToD(const std::array<float, 3>& V) { return {V[0], V[1], V[2]}; }
FVec3d DeltaRotate(const FBoneView& Off, const FBoneView& On, const FVec3d& V)
{
    FQuatd A{Off.Rotation[0], -Off.Rotation[1], -Off.Rotation[2], -Off.Rotation[3]};
    FQuatd B{On.Rotation[0], On.Rotation[1], On.Rotation[2], On.Rotation[3]};
    for(FQuatd* Q : {&A, &B})
    {
        double Length = 0.;
        for(double Value : *Q) Length += Value * Value;
        Length = std::sqrt(Length);
        for(double& Value : *Q) Value /= Length;
    }
    return RotateD(B, RotateD(A, V));
}
double FootFkError(const FFootStep& Off, const FFootStep& On, int Side)
{
    const auto& A = Off.Feet[Side];
    const auto& B = On.Feet[Side];
    const int Start = 1 + Side * 4;
    const std::array<FVec3d, 3> EdgesOff{SubD(ToD(A.KneePosition), ToD(A.HipPosition)), SubD(ToD(A.AnklePosition), ToD(A.KneePosition)),
                                       SubD(ToD(A.ToePosition), ToD(A.AnklePosition))};
    const std::array<FVec3d, 3> EdgesOn{SubD(ToD(B.KneePosition), ToD(B.HipPosition)), SubD(ToD(B.AnklePosition), ToD(B.KneePosition)),
                                      SubD(ToD(B.ToePosition), ToD(B.AnklePosition))};
    double Error = 0.;
    for(int Index = 0; Index < 3; ++Index)
        Error = std::max(Error, DistanceD(EdgesOn[Index], DeltaRotate(Off.Bones[Start + Index], On.Bones[Start + Index], EdgesOff[Index])));
    Error = std::max(Error, DistanceD(ToD(On.Bones[Start + 2].Position), ScaleD(AddD(ToD(B.AnklePosition), ToD(B.ToePosition)), 0.5)));
    Error = std::max(Error, DistanceD(SubD(ToD(On.Bones[Start + 3].Position), ToD(B.ToePosition)),
                                     DeltaRotate(Off.Bones[Start + 3], On.Bones[Start + 3], SubD(ToD(Off.Bones[Start + 3].Position), ToD(A.ToePosition)))));
    return Error;
}
void AnalyzeFootRuns(const std::string& Name, const std::vector<FFootStep>& Off, const std::vector<FFootStep>& On,
                     const std::vector<uint8_t>& Labels)
{
    bool bSameMotion = Off.size() == On.size();
    bool bDisplayLabel = true, bLabelBits = true, bOffIdentity = true, bOffUnlocked = true, bAnimShared = true, bLockLabel = true;
    bool bSearchConsistent = true, bFirstLocked = true;
    bool bUnchangedBody = true, bFreeIdentity = true, bFiniteFeet = true, bFixedLock = true, bToeRotation = true;
    int LiftedSteps = 0, BlendedLiftSteps = 0;
    double MaxFkError = 0., MaxNormErrorOn = 0., MaxNormErrorOff = 0.;
    int Searches = 0, Fallbacks = 0, Filtered = 0;
    float MaxBoneError = 0.f;
    double IkMsSum = 0., IkMsMax = 0.;
    for(size_t Index = 0; bSameMotion && Index < On.size(); ++Index)
    {
        const FFootStep& A = Off[Index];
        const FFootStep& B = On[Index];
        bSameMotion &= A.Root == B.Root && A.Yaw == B.Yaw && A.Source == B.Source && A.SearchCount == B.SearchCount
            && A.Query == B.Query && A.Result == B.Result && A.Contacts == B.Contacts && A.bFallback == B.bFallback
            && A.bBlending == B.bBlending;
        if(!B.bBlending) bDisplayLabel &= B.Source >= 0 && B.Source < static_cast<int>(Labels.size()) && B.Contacts == Labels[B.Source];
        for(int Bone = 0; Bone < SceneCapsuleCount; ++Bone)
        {
            if(Bone == 0 || Bone >= 9) bUnchangedBody &= SameBone(A.Bones[Bone], B.Bones[Bone]);
            double Norm = 0.;
            for(float V : B.Bones[Bone].Position) bFiniteFeet &= std::isfinite(V);
            for(float V : B.Bones[Bone].Rotation) Norm += static_cast<double>(V) * V;
            bFiniteFeet &= std::isfinite(Norm);
            // 存量动画/胶囊补偿不保证全部单位化；只对本次实际修改的三块骨检查单位误差。
            if(Bone >= 1 && Bone <= 8 && (Bone - 1) % 4 < 3)
            {
                const auto& Foot = B.Feet[(Bone - 1) / 4];
                if(Foot.bLocked || Foot.bFading) bFiniteFeet &= std::fabs(Norm - 1.) <= 1e-5;
            }
            MaxNormErrorOn = std::max(MaxNormErrorOn, std::fabs(Norm - 1.));
            double OffNorm = 0.;
            for(float V : A.Bones[Bone].Rotation) OffNorm += static_cast<double>(V) * V;
            MaxNormErrorOff = std::max(MaxNormErrorOff, std::fabs(OffNorm - 1.));
        }
        for(int Side = 0; Side < 2; ++Side)
        {
            const FFootView& F = A.Feet[Side];
            const FFootView& G = B.Feet[Side];
            const bool bBit = ((B.Contacts >> Side) & 1) != 0;
            bLabelBits &= F.bContactLabel == bBit && G.bContactLabel == bBit;
            bOffIdentity &= F.ToePosition == F.AnimToePosition;
            bOffUnlocked &= !F.bLocked && !F.bFading;
            bAnimShared &= G.AnimToePosition == F.ToePosition && G.HipPosition == F.HipPosition;
            bLockLabel &= !G.bLocked || G.bContactLabel;
            bFiniteFeet &= std::isfinite(G.ClampError) && G.ClampError >= 0.f && std::isfinite(G.LiftAngle)
                && G.LiftAngle >= 0.f && G.LiftAngle <= Pi / 4.f + 1e-6f && G.InvalidSteps == 0;
            for(const auto& Position : {G.HipPosition, G.KneePosition, G.AnklePosition, G.ToePosition, G.AnimToePosition})
                for(float Value : Position) bFiniteFeet &= std::isfinite(Value);
            MaxFkError = std::max(MaxFkError, FootFkError(A, B, Side));
            const int Start = 1 + Side * 4;
            bToeRotation &= A.Bones[Start + 3].Rotation == B.Bones[Start + 3].Rotation;
            if(!G.bLocked && !G.bFading)
                for(int Bone = Start; Bone < Start + 4; ++Bone) bFreeIdentity &= SameBone(A.Bones[Bone], B.Bones[Bone]);
            LiftedSteps += G.bLocked && G.LiftAngle > 1e-5f;
            BlendedLiftSteps += G.bLocked && G.LiftAngle > 1e-5f && B.bBlending;
            if(Index > 0 && G.bLocked && On[Index - 1].Feet[Side].bLocked)
                bFixedLock &= G.LockPosition == On[Index - 1].Feet[Side].LockPosition;
            if(Index == 0) bFirstLocked &= G.bLocked == G.bContactLabel || G.ReleaseCounts[static_cast<int>(EFootRelease::Clamp)] > 0;
            MaxBoneError = std::max({MaxBoneError,
                                     std::fabs(Distance3(G.HipPosition, G.KneePosition) - Distance3(F.HipPosition, F.KneePosition)),
                                     std::fabs(Distance3(G.KneePosition, G.AnklePosition) - Distance3(F.KneePosition, F.AnklePosition)),
                                     std::fabs(Distance3(G.AnklePosition, G.ToePosition) - Distance3(F.AnklePosition, F.ToePosition))});
        }
        if(Index > 0)
        {
            IkMsSum += B.IkMs;
            IkMsMax = std::max(IkMsMax, B.IkMs);
            if(B.SearchCount != On[Index - 1].SearchCount)
            {
                ++Searches;
                Fallbacks += B.bFallback;
                if(B.Query >= 0 && !B.bFallback) { ++Filtered; bSearchConsistent &= B.Result == B.Query; }
            }
        }
    }
    Require(bSameMotion, Name + ": IK on/off root, yaw, source frame and search sequence identical");
    const bool bFixed = Name == "foot fixed";
    Require(FootRunDigest(Off, true) == (bFixed ? 0xb89ab4eaa7ab55c2ull : 0x0421f8d7be9131e4ull)
            && FootRunDigest(On, false) == (bFixed ? 0x6f9818fbcc5da390ull : 0xd0268e226a0d47a6ull),
            Name + ": IK-off pose and motion sequence equal pre-pivot hash anchors");
    Require(bUnchangedBody && bFreeIdentity && bToeRotation, Name + ": pelvis/body untouched, free legs identical, toe retains animation rotation");
    std::cout << Name << " validation finite=" << bFiniteFeet << " fixedLock=" << bFixedLock << " fk=" << MaxFkError
              << " invalid=" << On.back().Feet[0].InvalidSteps << '/' << On.back().Feet[1].InvalidSteps
              << " normError off=" << MaxNormErrorOff << " on=" << MaxNormErrorOn << '\n';
    Require(bFiniteFeet && bFixedLock && MaxFkError <= 1e-4, Name + ": finite unit rotations, fixed lock points, capsule/chain FK <=1e-4m");
    Require(LiftedSteps > 0 && BlendedLiftSteps > 0, Name + ": actual heel lift exercised in locked and blended poses");
    auto Corrupt = On.front();
    Corrupt.Feet[0].ToePosition[0] += 0.002f;
    Require(FootFkError(Off.front(), Corrupt, 0) > 1e-4, Name + ": negative control displaced toe fails FK check");
    Corrupt = Off.front();
    Corrupt.Bones[1].Position[0] = std::nextafter(Corrupt.Bones[1].Position[0], std::numeric_limits<float>::infinity());
    Require(FootRunDigest({Corrupt}, true) != FootRunDigest({Off.front()}, true), Name + ": negative control float-bit change alters pose digest");
    Require(bDisplayLabel && bLabelBits, Name + ": displayed contacts equal labels[source] when not blending, per-foot bits match");
    Require(bSearchConsistent, Name + ": filtered searches return candidates with the query contact label");
    Require(bOffIdentity && bOffUnlocked, Name + ": IK off never locks and shows the animation toe exactly");
    Require(bAnimShared, Name + ": IK on keeps hip and animation toe equal to the IK-off pose");
    Require(bLockLabel && bFirstLocked, Name + ": locks only during contact labels; reset locks feet already in contact");
    std::cout << Name << " searches=" << Searches << " filtered=" << Filtered << " fallbacks=" << Fallbacks
              << " maxBoneLengthChange=" << MaxBoneError << '\n';
    Require(MaxBoneError <= 1e-4f, Name + ": IK keeps bone lengths (float pose, <=1e-4m)");

    const FFootRunMetrics MetricsOff = MeasureFootRun(Off);
    const FFootRunMetrics MetricsOn = MeasureFootRun(On);
    const auto& Last = On.back().Feet;
    std::cout << Name << " offPoseHash=" << std::hex << FootRunDigest(Off, true) << " motionHash=" << FootRunDigest(On, false) << std::dec << '\n';
    for(int Side = 0; Side < 2; ++Side)
    {
        int ContactSteps = 0, LockedContactSteps = 0, ConsecutiveFade = 0, MaxConsecutiveFade = 0;
        bool bHandbacksHealthy = true;
        for(const auto& Step : On)
        {
            const auto& Foot = Step.Feet[Side];
            ContactSteps += Foot.bContactLabel;
            LockedContactSteps += Foot.bContactLabel && Foot.bLocked;
            ConsecutiveFade = Foot.bFading ? ConsecutiveFade + 1 : 0;
            MaxConsecutiveFade = std::max(MaxConsecutiveFade, ConsecutiveFade);
            bHandbacksHealthy &= Foot.ForcedHandbacks == 0 && Foot.MaxSettleSteps <= 30;
        }
        std::cout << Name << " side=" << Side << " coverage=" << LockedContactSteps << '/' << ContactSteps << " conflicts=" << Last[Side].ConflictSteps
                  << " fadeConflicts=" << Last[Side].FadeConflictSteps << " maxFkError=" << MaxFkError << " liftedSteps=" << LiftedSteps
                  << " maxConsecutiveFade=" << MaxConsecutiveFade << " maxSettle=" << Last[Side].MaxSettleSteps
                  << " forcedHandbacks=" << Last[Side].ForcedHandbacks << '\n';
        Require(bHandbacksHealthy && MaxConsecutiveFade <= 120,
                Name + ": no forced handback, settle <=30 steps and consecutive fade <=120 steps, side=" + std::to_string(Side));
        const int BaselineLocked = bFixed ? (Side == 0 ? 334 : 356) : (Side == 0 ? 13427 : 13362);
        const int BaselineEarly = bFixed ? (Side == 0 ? 2 : 1) : (Side == 0 ? 135 : 147);
        Budget(LockedContactSteps >= BaselineLocked && Last[Side].ReleaseCounts[2] + Last[Side].ReleaseCounts[3] <= BaselineEarly,
               Name + " lock coverage/early releases not degraded, side=" + std::to_string(Side),
               "locked=" + std::to_string(LockedContactSteps) + " baseline=" + std::to_string(BaselineLocked));
    }
    std::cout << Name << " slip off=" << MetricsOff.LabelSlip << "m on=" << MetricsOn.LabelSlip << "m ratio="
              << (MetricsOff.LabelSlip > 0. ? MetricsOn.LabelSlip / MetricsOff.LabelSlip : 0.)
              << " lockedSteps=" << MetricsOn.LockedSteps << " maxLockedStep=" << MetricsOn.MaxLockedStep
              << " maxLockOffset=" << MetricsOn.MaxLockOffset << " maxLockHeightError=" << MetricsOn.MaxLockHeightError
              << " kneeDelta off=" << MetricsOff.MaxKneeDelta << " on=" << MetricsOn.MaxKneeDelta
              << " extra max=" << MetricsOn.MaxExtra << " fadeSteps=" << MetricsOn.FadeSteps << " fadeExtra=" << MetricsOn.MaxFadeExtra
              << " ikMs mean=" << (On.size() > 1 ? IkMsSum / (On.size() - 1) : 0.) << " max=" << IkMsMax << '\n';
    std::cout << Name << " diag slipLocked=" << MetricsOn.LabelSlipLocked << " lockedOverClamped=" << MetricsOn.LockedOverClamped << " maxLockedStepUnclamped=" << MetricsOn.MaxLockedStepUnclamped
              << " kneeByState locked/fade/switch/free=" << MetricsOn.MaxKneeByState[0] << '/' << MetricsOn.MaxKneeByState[1] << '/'
              << MetricsOn.MaxKneeByState[2] << '/' << MetricsOn.MaxKneeByState[3] << " clampAtMax=" << MetricsOn.MaxClampByKnee[0] << '/'
              << MetricsOn.MaxClampByKnee[1] << '/' << MetricsOn.MaxClampByKnee[2] << '/' << MetricsOn.MaxClampByKnee[3] << '\n';
    for(int Side = 0; Side < 2; ++Side)
        std::cout << Name << (Side == 0 ? " left" : " right") << " releases label=" << Last[Side].ReleaseCounts[1] << " distance="
                  << Last[Side].ReleaseCounts[2] << " clamp=" << Last[Side].ReleaseCounts[3] << " disabled=" << Last[Side].ReleaseCounts[4] << '\n';
    Budget(MetricsOn.LockedSteps > 0 && MetricsOn.LockedOver == 0, Name + " locked toe horizontal step <=1e-3m (100%)",
           "steps=" + std::to_string(MetricsOn.LockedSteps) + " over=" + std::to_string(MetricsOn.LockedOver)
           + " max=" + std::to_string(MetricsOn.MaxLockedStep));
    if(MetricsOff.LabelSlip >= 1.)
        Budget(MetricsOn.LabelSlip <= 0.2 * MetricsOff.LabelSlip, Name + " contact-label slip IK on <=20% of IK off",
               "off=" + std::to_string(MetricsOff.LabelSlip) + " on=" + std::to_string(MetricsOn.LabelSlip));
    else
        std::cout << "NOT JUDGED: " << Name << " slip ratio (IK off slip " << MetricsOff.LabelSlip << "m < 1m)\n";
    // 保留最坏释放步的归因证据，不仅给出总数。
    std::vector<std::pair<float, std::pair<size_t, int>>> Outliers;
    for(size_t Index = 1; Index < On.size(); ++Index)
        for(int Side = 0; Side < 2; ++Side)
        {
            const auto& P = On[Index - 1].Feet[Side];
            const auto& C = On[Index].Feet[Side];
            const float Extra = Distance3(C.ToePosition, P.ToePosition) - Distance3(C.AnimToePosition, P.AnimToePosition);
            if(Extra > 0.02f) Outliers.push_back({Extra, {Index, Side}});
        }
    std::sort(Outliers.rbegin(), Outliers.rend());
    for(size_t Rank = 0; Rank < std::min(size_t(3), Outliers.size()); ++Rank)
    {
        const size_t Index = Outliers[Rank].second.first;
        const int Side = Outliers[Rank].second.second;
        const auto& P = On[Index - 1].Feet[Side];
        const auto& C = On[Index].Feet[Side];
        std::cout << Name << " outlier step=" << Index << " side=" << Side << " extra=" << Outliers[Rank].first
                  << " lift=" << Degrees(P.LiftAngle) << "->" << Degrees(C.LiftAngle) << " residual=" << P.ClampError << "->" << C.ClampError
                  << " knee=" << KneeDegrees(P) << "->" << KneeDegrees(C) << " lock=" << P.bLocked << C.bLocked << " fade=" << P.bFading << C.bFading
                  << " toeGap=" << Distance3(P.ToePosition, P.AnimToePosition) << "->" << Distance3(C.ToePosition, C.AnimToePosition)
                  << " hipStep=" << Distance3(P.HipPosition, C.HipPosition) << '\n';
    }
    const float KneeLimit = std::max(15.f, MetricsOff.MaxKneeDelta);
    Budget(MetricsOn.MaxKneeDelta <= KneeLimit, Name + " knee angle step <=max(15deg, IK off max)",
           "on=" + std::to_string(MetricsOn.MaxKneeDelta) + " off=" + std::to_string(MetricsOff.MaxKneeDelta));
    Budget(MetricsOn.MaxKneeDelta <= 15.f, Name + " original absolute knee angle step <=15deg", "on=" + std::to_string(MetricsOn.MaxKneeDelta));
    Budget(MetricsOn.ExtraOver == 0, Name + " display toe step <= animation toe step + 0.02m",
           "over=" + std::to_string(MetricsOn.ExtraOver) + " max=" + std::to_string(MetricsOn.MaxExtra)
           + " fadeMax=" + std::to_string(MetricsOn.MaxFadeExtra));
}

bool AnyLocked(const FSceneSnapshot& Snapshot) { return Snapshot.Feet[0].bLocked || Snapshot.Feet[1].bLocked; }
bool SameFeet(const FSceneSnapshot& A, const FSceneSnapshot& B)
{
    bool bSame = true;
    for(int Side = 0; Side < 2; ++Side)
    {
        const FFootView& F = A.Feet[Side];
        const FFootView& G = B.Feet[Side];
        bSame &= F.ToePosition == G.ToePosition && F.KneePosition == G.KneePosition && F.bLocked == G.bLocked
            && F.bFading == G.bFading && F.LockPosition == G.LockPosition && F.ReleaseCounts == G.ReleaseCounts
            && F.MaxSettleSteps == G.MaxSettleSteps && F.ForcedHandbacks == G.ForcedHandbacks;
        for(int Bone = 1 + Side * 4; Bone < 5 + Side * 4; ++Bone) bSame &= SameBone(A.Bones[Bone], B.Bones[Bone]);
    }
    return bSame;
}

void ProbeFootIkSettings(FSandboxScene& Scene)
{
    FSceneCommand Restore;
    Restore.bRestoreDefaults = true;
    Restore.bReset = true;
    SetMode(Scene, ESceneMode::MotionMatching);
    Scene.ApplyCommand(Restore);
    SetPaused(Scene, false);
    Require(Scene.GetSettings().bFootIk, "foot IK enabled by default");
    const auto First = Scene.GetSnapshot();
    bool bFirst = true;
    for(const FFootView& Foot : First.Feet) bFirst &= Foot.bLocked == Foot.bContactLabel || Foot.LastRelease == EFootRelease::Clamp;
    Require(bFirst, "matching reset locks feet already in contact");

    // 暂停中关闭IK：快照不变，下一步已锁定的脚以Disabled释放并淡出；淡出结束后显示动画姿态且不再锁定。
    for(int Step = 0; Step < 240; ++Step)
    {
        Scene.Tick({0.f, -1.f, 0.f, 0.f});
        if(Step >= 30 && AnyLocked(Scene.GetSnapshot())) break;
    }
    const auto Locked = Scene.GetSnapshot();
    Require(AnyLocked(Locked), "walking locks a contact foot");
    SetPaused(Scene, true);
    FSceneSettings Off = Scene.GetSettings();
    Off.bFootIk = false;
    Require(Scene.SetSettings(Off), "foot IK disable accepted while paused");
    const auto Frozen = Scene.GetSnapshot();
    Require(SameFeet(Frozen, Locked) && Frozen.RootPosition == Locked.RootPosition && Frozen.NextFrame == Locked.NextFrame,
            "paused IK toggle changes nothing until the next step");
    FSceneCommand Single;
    Single.bSingleStep = true;
    Scene.ApplyCommand(Single);
    Scene.Tick({0.f, -1.f, 0.f, 0.f});
    const auto Released = Scene.GetSnapshot();
    bool bReleased = true;
    const int Disabled = static_cast<int>(EFootRelease::Disabled);
    for(int Side = 0; Side < 2; ++Side)
        if(Locked.Feet[Side].bLocked)
            bReleased &= !Released.Feet[Side].bLocked && Released.Feet[Side].bFading && Released.Feet[Side].LastRelease == EFootRelease::Disabled
                && Released.Feet[Side].ReleaseCounts[Disabled] == Locked.Feet[Side].ReleaseCounts[Disabled] + 1;
    Require(bReleased && Released.NextFrame == Locked.NextFrame + 1, "next step releases locked feet as Disabled and fades");
    SetPaused(Scene, false);
    bool bNeverLocked = true, bAnimated = true;
    int FadeEnd = -1;
    for(int Step = 0; Step < 120; ++Step)
    {
        Scene.Tick({0.f, -1.f, 0.f, 0.f});
        const auto Now = Scene.GetSnapshot();
        bNeverLocked &= !AnyLocked(Now);
        if(FadeEnd < 0 && !Now.Feet[0].bFading && !Now.Feet[1].bFading) FadeEnd = Step;
        if(FadeEnd >= 0)
            for(const FFootView& Foot : Now.Feet) bAnimated &= !Foot.bFading && Foot.ToePosition == Foot.AnimToePosition;
    }
    std::cout << "foot IK disable fade ended " << FadeEnd + 1 << " steps after the release step\n";
    Require(bNeverLocked && FadeEnd >= 0 && bAnimated, "IK off never locks; after fade the display equals the animation");
    for(const auto& Foot : Scene.GetSnapshot().Feet)
        Require(Foot.ForcedHandbacks == 0 && Foot.MaxSettleSteps <= 30, "near-origin disable completes without forced handback");

    // Reset与模式切换保留关闭状态并清空锁定；RestoreDefaults恢复开启；模式切换即重置且可复现。
    FSceneCommand Reset;
    Reset.bReset = true;
    Scene.ApplyCommand(Reset);
    const auto AfterReset = Scene.GetSnapshot();
    bool bCleared = !Scene.GetSettings().bFootIk && AfterReset.NextFrame == 0;
    for(const FFootView& Foot : AfterReset.Feet)
        bCleared &= !Foot.bLocked && !Foot.bFading && Foot.ReleaseCounts == std::array<int, FootReleaseKinds>{}
            && Foot.MaxSettleSteps == 0 && Foot.ForcedHandbacks == 0 && Foot.ToePosition == Foot.AnimToePosition;
    Require(bCleared, "reset keeps IK off and clears lock state and release counts");
    SetMode(Scene, ESceneMode::Physics);
    SetMode(Scene, ESceneMode::MotionMatching);
    Require(!Scene.GetSettings().bFootIk, "mode switch keeps IK off");
    Scene.ApplyCommand(Restore);
    Require(Scene.GetSettings().bFootIk && SameFeet(Scene.GetSnapshot(), First), "restore defaults re-enables foot IK and reproduces reset pose");
    SetPaused(Scene, false);
    for(int Step = 0; Step < 90; ++Step) Scene.Tick({1.f, 0.f, 0.f, 0.f});
    SetMode(Scene, ESceneMode::Physics);
    SetMode(Scene, ESceneMode::MotionMatching);
    const auto Switched = Scene.GetSnapshot();
    Require(Switched.NextFrame == 0 && SameFeet(Switched, First), "mode switch resets matching and foot locks reproducibly");
}

// 用真实MM行走把世界坐标推进到100m，覆盖float舍入后关闭IK的交还边界；不修改root或检索状态。
void ProbeFootIkFarDisable(FSandboxScene& Scene)
{
    FSceneCommand Restore;
    Restore.bRestoreDefaults = true;
    Restore.bReset = true;
    SetMode(Scene, ESceneMode::MotionMatching);
    Require(Scene.ApplyCommand(Restore), "far-origin restore accepted");
    SetPaused(Scene, false);
    const auto Origin = Scene.GetSnapshot().RootPosition;
    FSceneSnapshot Locked;
    int TravelSteps = 0;
    for(; TravelSteps < 12000; ++TravelSteps)
    {
        Scene.Tick({0.f, -1.f, 0.f, 0.f});
        Locked = Scene.GetSnapshot();
        if(HorizontalDistance(Locked.RootPosition, Origin) >= 100.f && AnyLocked(Locked)) break;
    }
    Require(TravelSteps < 12000 && HorizontalDistance(Locked.RootPosition, Origin) >= 100.f && AnyLocked(Locked),
            "far-origin disable reaches >=100m with a locked foot");
    auto Settings = Scene.GetSettings();
    Settings.bFootIk = false;
    Require(Scene.SetSettings(Settings), "far-origin IK disable accepted");
    int FadeEnd = -1;
    bool bAnimated = true, bNeverLocked = true, bReleased = true, bHandbacksHealthy = true;
    std::array<int, 2> ConsecutiveFade{}, MaxConsecutiveFade{};
    const int Disabled = static_cast<int>(EFootRelease::Disabled);
    for(int Step = 0; Step < 120; ++Step)
    {
        Scene.Tick({0.f, -1.f, 0.f, 0.f});
        const auto Now = Scene.GetSnapshot();
        bNeverLocked &= !AnyLocked(Now);
        if(FadeEnd < 0 && !Now.Feet[0].bFading && !Now.Feet[1].bFading) FadeEnd = Step;
        for(int Side = 0; Side < 2; ++Side)
        {
            const auto& Foot = Now.Feet[Side];
            if(Step == 0 && Locked.Feet[Side].bLocked)
                bReleased &= Foot.bFading && Foot.LastRelease == EFootRelease::Disabled
                    && Foot.ReleaseCounts[Disabled] == Locked.Feet[Side].ReleaseCounts[Disabled] + 1;
            ConsecutiveFade[Side] = Foot.bFading ? ConsecutiveFade[Side] + 1 : 0;
            MaxConsecutiveFade[Side] = std::max(MaxConsecutiveFade[Side], ConsecutiveFade[Side]);
            bHandbacksHealthy &= Foot.ForcedHandbacks == 0 && Foot.MaxSettleSteps <= 30 && Foot.InvalidSteps == 0;
            if(FadeEnd >= 0) bAnimated &= !Foot.bFading && Foot.ToePosition == Foot.AnimToePosition;
        }
    }
    const auto Last = Scene.GetSnapshot();
    std::cout << "foot far disable distance=" << HorizontalDistance(Locked.RootPosition, Origin)
              << " travelSteps=" << TravelSteps + 1 << " fadeEnd=" << FadeEnd + 1
              << " maxConsecutiveFade=" << MaxConsecutiveFade[0] << '/' << MaxConsecutiveFade[1]
              << " maxSettle=" << Last.Feet[0].MaxSettleSteps << '/' << Last.Feet[1].MaxSettleSteps
              << " forcedHandbacks=" << Last.Feet[0].ForcedHandbacks << '/' << Last.Feet[1].ForcedHandbacks << '\n';
    Require(bReleased && bNeverLocked && FadeEnd >= 0 && bAnimated && bHandbacksHealthy,
            "far-origin disable finishes within 120 steps, toes exactly animation, no forced handback or invalid solve");
    auto Corrupt = Last;
    Corrupt.Feet[0].ToePosition[0] = std::nextafter(Corrupt.Feet[0].AnimToePosition[0], std::numeric_limits<float>::infinity());
    Require(Corrupt.Feet[0].ToePosition != Corrupt.Feet[0].AnimToePosition,
            "negative control one-ULP displaced toe fails exact animation handback check");
}

void ProbeFootIk(FSandboxScene& Scene)
{
    uint64_t SourceHash = 0;
    std::string Reason;
    FFootContactFile File;
    Require(HashFileFnv1a64(GetSandboxResourcePath("long_motion.bin"), SourceHash, Reason)
            && LoadFootContacts(GetSandboxResourcePath("foot_contacts.bin"), Seams[6], SourceHash, File, Reason) == EFootContactError::None,
            "probe loads foot contact labels");
    ProbeFootIkSettings(Scene);
    ProbeFootIkFarDisable(Scene);
    FSceneCommand Restore;
    Restore.bRestoreDefaults = true;
    Restore.bReset = true;
    SetMode(Scene, ESceneMode::MotionMatching);
    Scene.ApplyCommand(Restore);
    SetPaused(Scene, false);
    struct FScript
    {
        const char* Name;
        int Steps;
        std::function<FInputScript()> Make;
    };
    for(const FScript& Script : {FScript{"foot fixed", 630, [] { return FixedFootScript(); }},
                                 FScript{"foot random", 30000, [] { return RandomFootScript(20240531u); }}})
    {
        const auto Off = RunFootScript(Scene, false, Script.Steps, Script.Make());
        const auto On = RunFootScript(Scene, true, Script.Steps, Script.Make());
        AnalyzeFootRuns(Script.Name, Off, On, File.Labels);
        const auto Replay = RunFootScript(Scene, true, Script.Steps, Script.Make());
        Require(FootRunDigest(On, true) == FootRunDigest(Replay, true), std::string(Script.Name) + ": reset/replay pose digest identical");
    }
    Scene.ApplyCommand(Restore);
    SetPaused(Scene, false);
}
}

int main(int ArgCount, char** Arguments)
{
    const std::string Option = ArgCount == 3 ? Arguments[2] : "";
    if(ArgCount != 2 && !(ArgCount == 3 && (Option == "--digest-only" || Option == "--ground-perf" || Option == "--foot")))
    { std::cerr << "usage: SceneProbe OUTPUT_PREFIX [--digest-only|--ground-perf|--foot]\n"; return 2; }
    try
    {
        if(Option == "--digest-only")
        {
            FSandboxScene Scene;
            PrintDigests(Scene);
            return 0;
        }
        // 工具输出路径在场景创建前解析为绝对路径。
        const std::string Prefix = std::filesystem::absolute(Arguments[1]).string();
        if(Option == "--foot")
        {
            // 只跑无GPU的标签/IK/滑步探针；完整运行同样包含这些检查。
            FSandboxScene Scene;
            ProbeTwoBoneIk();
            ProbeFootPivotIk();
            ProbeFootContactFile(Scene, Prefix);
            ProbeFootIk(Scene);
            std::cout << "foot hard checks=PASS budgetFailures=" << BudgetFailures << '\n';
            return BudgetFailures ? 3 : 0;
        }
        FAppWindow Window({720, 720, "SceneProbe", true});
        FSandboxScene Scene;
        FGpuRenderer Gpu(Window, 720, 720);
        Scene.InitializeGpu(Gpu);
        if(Option == "--ground-perf")
        {
            // 仅工具中的同姿态开/关对照：排除新旧BRDF与MM差异，不新增产品CLI。
            // 光源模型始终保留；计时只包含离屏绘制，不含呈现/交换，GPU timestamp异步滞后。
            for(int Frame = 0; Frame <= 60; ++Frame) Scene.Tick({});
            auto Settings = Scene.GetSettings();
            Settings.bShowCharacter = Settings.bShowCapsules = Settings.bShowContacts = false;
            SetPaused(Scene, true);
            for(bool bCpu : {false, true})
            {
                for(int Run = 1; Run <= 3; ++Run)
                {
                    for(bool bVisible : {Run % 2 == 0, Run % 2 != 0})
                    {
                        Settings.Ground.bVisible = bVisible;
                        Require(Scene.SetSettings(Settings), "ground performance toggle accepted");
                        for(int Frame = 0; Frame < 600; ++Frame)
                        {
                            Window.BeginFrame();
                            const auto Start = std::chrono::steady_clock::now();
                            if(bCpu) Scene.RenderCpu();
                            else Scene.BuildRenderData(Gpu);
                            const double Cost = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - Start).count();
                            Window.EndFrame();
                            std::cout << "GROUND_PERF " << (bCpu ? "cpu" : "gpu") << ' ' << Run << ' ' << bVisible << ' ' << Frame
                                      << ' ' << Cost << ' ' << (bCpu ? 0. : Gpu.GetLastGpuMs()) << '\n';
                        }
                    }
                }
            }
            return 0;
        }
        for(int Frame = 0; Frame <= 60; ++Frame) Scene.Tick({});
        Scene.BuildRenderData(Gpu);
        FImageRGB8 Reference;
        Gpu.ReadbackRGB8(Reference);
        Save(Reference, Prefix + "_gpu.png");
        Scene.RenderCpu();
        FImageRGB8 Cpu;
        Scene.ReadCpuRGB8(Cpu);
        Save(Cpu, Prefix + "_cpu.png");
        FSceneCommand Pause;
        Pause.bTogglePause = true;
        Scene.ApplyCommand(Pause);
        const auto Before = Scene.GetSnapshot();
        for(int Frame = 0; Frame < 20; ++Frame) Scene.Tick({1.f, -1.f, 1.f, 1.f});
        const auto After = Scene.GetSnapshot();
        Require(After.NextFrame == Before.NextFrame && After.RootPosition == Before.RootPosition
                && After.LightPosition == Before.LightPosition && After.LightRadiance == Before.LightRadiance, "pause freezes state");
        FSceneCommand Step;
        Step.bSingleStep = true;
        Scene.ApplyCommand(Step);
        Scene.Tick({});
        Scene.Tick({});
        Require(Scene.GetSnapshot().NextFrame == Before.NextFrame + 1 && Scene.GetSnapshot().bPaused, "single step exactly once");
        auto Settings = Scene.GetSettings();
        Settings.Gravity = std::numeric_limits<double>::quiet_NaN();
        Require(!Scene.SetSettings(Settings) && Scene.GetSettings().Gravity == 9.8, "reject nonfinite settings");
        Settings = Scene.GetSettings();
        Settings.MatchingWeights = {0.f, 0.f, 0.f};
        Require(!Scene.SetSettings(Settings), "reject all zero weights");
        Settings.MatchingWeights = {1.f, -1.f, 1.f};
        Require(!Scene.SetSettings(Settings), "reject negative weight");
        FSceneCommand Reset;
        Reset.bRestoreDefaults = true;
        Reset.bReset = true;
        Scene.ApplyCommand(Reset);
        Scene.ApplyCommand(Pause);
        for(int Frame = 0; Frame <= 60; ++Frame) Scene.Tick({});
        Scene.BuildRenderData(Gpu);
        FImageRGB8 Repeated;
        Gpu.ReadbackRGB8(Repeated);
        Require(Reference.Pixels == Repeated.Pixels, "reset reproduces GPU frame60");
        Scene.RenderCpu();
        Scene.ReadCpuRGB8(Repeated);
        Require(Cpu.Pixels == Repeated.Pixels, "reset reproduces CPU frame60");

        ProbeGroundGeometry();
        ProbeMaterialsAndGround(Scene, Gpu);
        ProbeRouter();
        ProbeMatching(Scene);
        ProbeTwoBoneIk();
        ProbeFootPivotIk();
        ProbeFootContactFile(Scene, Prefix);
        ProbeFootIk(Scene);
        const auto First = RunRandom(Scene, 30000);
        const auto Second = RunRandom(Scene, 30000);
        Require(First.Trace == Second.Trace, "random 30000 steps root and source deterministic");
        const auto Different = RunRandom(Scene, 1000, 20240532u);
        Require(!std::equal(Different.Trace.begin(), Different.Trace.end(), First.Trace.begin()),
                "negative control different seed changes trajectory");
        std::cout << "random searches=" << First.Searches << " avgSearchMs=" << (First.Searches ? First.SearchMs / First.Searches : 0.)
                  << " avgTickMs=" << First.TickMs << '\n';
        Budget(First.MaxDisplacement <= 0.1f, "step displacement <=0.1m", "max=" + std::to_string(First.MaxDisplacement));
        Budget(First.MaxYawDelta <= 20.f, "step yaw delta <=20deg", "max=" + std::to_string(First.MaxYawDelta) + " " + First.MaxYawWhere);

        // Physics -> MM -> Physics 后重新执行0..60，必须恢复历史帧60。
        FSceneCommand Back;
        Back.bSetMode = true;
        Back.Mode = ESceneMode::Physics;
        Require(Scene.ApplyCommand(Back) && !Scene.GetSnapshot().bPaused, "return to physics");
        for(int Frame = 0; Frame <= 60; ++Frame) Scene.Tick({});
        Scene.BuildRenderData(Gpu);
        Gpu.ReadbackRGB8(Repeated);
        Require(Reference.Pixels == Repeated.Pixels, "mode round trip reproduces GPU frame60");
        Scene.RenderCpu();
        Scene.ReadCpuRGB8(Repeated);
        Require(Cpu.Pixels == Repeated.Pixels, "mode round trip reproduces CPU frame60");
        ProbeFrameCaptures(Scene, Gpu, Prefix);
        std::cout << "NOT RUN: no-candidate injection, clip exhaustion, hidden device wiring, real UI\n";
        std::cout << "hard checks=PASS budgetFailures=" << BudgetFailures << '\n';
        return BudgetFailures ? 3 : 0;
    }
    catch(...) { return 1; }
}
