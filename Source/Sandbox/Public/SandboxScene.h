#pragma once
#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include "SceneTypes.h"
#include "RenderTypes.h"

class FGpuRenderer;
// 旧对象只存在于PImpl，GL资源由调用方持有且先于窗口销毁。
class FSandboxScene
{
public:
    // FootContactsPath为空时读资源目录的foot_contacts.bin；标签缺失或校验失败时MM不可用、Physics照常。
    explicit FSandboxScene(const std::string& FootContactsPath = std::string());
    ~FSandboxScene();
    FSandboxScene(const FSandboxScene&) = delete;
    FSandboxScene& operator=(const FSandboxScene&) = delete;
    const FSceneSettings& GetSettings() const;
    bool SetSettings(const FSceneSettings& Settings);
    bool ApplyCommand(const FSceneCommand& Command);
    void Reset();
    void Tick(const FSceneInput& Input);
    FSceneSnapshot GetSnapshot() const;
    void InitializeGpu(FGpuRenderer& Renderer);
    void BuildRenderData(FGpuRenderer& Renderer);
    FImageRGBFView RenderCpu();
    void ReadCpuRGB8(FImageRGB8& Image);
private:
    struct FImpl;
    std::unique_ptr<FImpl> Impl;
};

// 资源目录下的文件路径（与场景加载同一定位规则）。
std::string GetSandboxResourcePath(const char* Name);
// 逐脚触地标签的生成统计；比例均为帧占比。
struct FFootContactStats
{
    int Frames = 0;
    std::array<double, 2> ContactRate = {};
    double BothRate = 0.;
    double FlightRate = 0.;
    // 镜像片段对(0,1)(2,3)(4,5)上“片段A左=片段B右且A右=B左”的帧占比
    std::array<double, 3> MirrorAgreement = {};
    // root水平速度 < 0.1m/s 的帧数及其双脚着地率
    int StillFrames = 0;
    double StillBothRate = 0.;
    // 与database单支撑脚标签（0左、1右）的一致率：支撑脚在新标签中为触地
    double SupportAgreement = 0.;
    int MinRun = 0;
    int MedianRun = 0;
};
// 由long_motion.bin前向运动学重算逐脚触地标签并序列化为foot_contacts.bin字节；Report为可读统计或失败原因。
bool BuildFootContactsFile(std::vector<uint8_t>& OutBytes, FFootContactStats& Stats, std::string& Report);
// 同上并写出OutputPath。
bool BakeFootContacts(const std::string& OutputPath, std::string& Report);
