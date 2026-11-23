#pragma once
#include <array>
#include <cstdint>
#include <string>

// Editor只交换这些值；不引用旧骨架/矩阵类型或GL。
enum class ESceneMode : uint8_t { Physics, MotionMatching };
enum class EControlMode : uint8_t { Global, Local, Off };
inline constexpr int SceneCapsuleCount = 22;
// 标准单次散射 metallic-roughness；Roughness 为感知粗糙度（alpha = r^2）。
struct FSurfaceMaterial
{
    float Metallic = 0.f;
    float Roughness = 0.1f;
};
// 世界固定的程序化棋盘地面。Physics 下平面取 ContactPlane；MM 下为 y = MatchingHeight 的水平面。
struct FGroundSettings
{
    bool bVisible = true;
    float MatchingHeight = 0.f;
    float CellSize = 0.5f;
    FSurfaceMaterial Material{0.f, 0.8f};
};
struct FSceneSettings
{
    double Gravity = 9.8;
    double Damping = 0.25;
    std::array<double, 4> ContactPlane = {0., 1., 0., 0.05};
    double ContactEpsilon = 1e-2;
    int Substeps = 1;
    EControlMode Control = EControlMode::Global;
    int ControlEndFrame = 900;
    float CameraDistance = 2.5f;
    float CameraHeight = 0.35f;
    float CameraSensitivity = 0.1f;
    bool bFollowLight = true;
    std::array<float, 3> LightPosition = {0.f, 1.f, 1.5f};
    float LightBase = 20.f;
    float LightAmplitude = 10.f;
    FSurfaceMaterial CharacterMaterial;
    std::array<FSurfaceMaterial, SceneCapsuleCount> CapsuleMaterials{};
    FGroundSettings Ground;
    bool bShowCapsules = true;
    bool bShowContacts = true;
    bool bShowCharacter = false;
    std::array<int, 3> LegacyWeights = {8, 3, 5};
    float MoveSpeed = 3.f;
    std::array<float, 3> MatchingWeights = {1.f, 1.f, 0.5f};
    // 脚锁IK：仅Motion Matching生效；关闭时已锁定的脚按释放流程淡出。
    bool bFootIk = true;
};
struct FSceneInput
{
    float MoveX = 0.f;
    float MoveY = 0.f;
    float LookX = 0.f;
    float LookY = 0.f;
};
struct FSceneCommand
{
    bool bSetMode = false;
    ESceneMode Mode = ESceneMode::Physics;
    bool bTogglePause = false;
    bool bSingleStep = false;
    bool bReset = false;
    bool bRestoreDefaults = false;
};
struct FBoneView
{
    std::array<float, 3> Position = {};
    std::array<float, 4> Rotation = {};
    bool bContact = false;
};
// 脚锁释放原因：标签下降沿、锁点远离动画脚尖、目标不可达、IK关闭。
enum class EFootRelease : uint8_t { None, Label, Distance, Clamp, Disabled };
inline constexpr int FootReleaseKinds = 5;
// 单脚（序号0左、1右）的只读状态；位置均为世界系，显示值已含IK。
struct FFootView
{
    std::array<float, 3> ToePosition = {};
    std::array<float, 3> AnimToePosition = {};
    std::array<float, 3> LockPosition = {};
    std::array<float, 3> HipPosition = {};
    std::array<float, 3> KneePosition = {};
    std::array<float, 3> AnklePosition = {};
    bool bContactLabel = false;
    bool bLocked = false;
    bool bFading = false;
    // 本步脚尖目标的实际FK残差（m），包含膝角优先的几何冲突及数值舍入；是否冲突看bConflict。
    float ClampError = 0.f;
    // 抬跟角（弧度）及膝角优先导致的几何冲突统计；仅供诊断，不改变锁定判定。
    float LiftAngle = 0.f;
    bool bConflict = false;
    int ConflictSteps = 0;
    int FadeConflictSteps = 0;
    int InvalidSteps = 0;
    // 权重归零后最长等待步数与超限强制交还次数；正常运行后者应为0。
    int MaxSettleSteps = 0;
    int ForcedHandbacks = 0;
    EFootRelease LastRelease = EFootRelease::None;
    std::array<int, FootReleaseKinds> ReleaseCounts = {};
};
struct FSceneSnapshot
{
    ESceneMode Mode = ESceneMode::Physics;
    bool bPaused = false;
    int NextFrame = 0;
    int LastFrame = -1;
    bool bMatchingAvailable = false;
    std::string Message;
    std::array<float, 3> RootPosition = {};
    float RootYaw = 0.f;
    float RootSpeed = 0.f;
    std::array<float, 3> CameraPosition = {};
    float CameraYaw = 0.f;
    std::array<float, 3> LightPosition = {};
    std::array<float, 3> LightRadiance = {};
    std::array<FBoneView, SceneCapsuleCount> Bones;
    // 当前模式实际使用的地面平面 A*x+B*y+C*z+D=0（未归一化原始系数）
    std::array<double, 4> GroundPlane = {};
    int SourceFrame = -1;
    int ClipStart = -1;
    bool bBlending = false;
    int SearchFrame = -1;
    int SearchCount = 0;
    int CandidateCount = 0;
    int ExcludedCount = 0;
    float SearchCost = 0.f;
    double SearchMs = 0.;
    // 逐脚触地标签（bit0左、bit1右）：当前显示、最近一次检索的查询与结果；回退表示无同标签候选而未过滤。
    int FootContacts = 0;
    int SearchQueryContacts = -1;
    int SearchResultContacts = -1;
    bool bSearchContactFallback = false;
    std::array<FFootView, 2> Feet;
    double FootIkMs = 0.;
};
struct FFrameStats
{
    uint64_t AppFrame = 0;
    double CostMs = 0.;
    double GpuMs = 0.;
    double IntervalMs = 0.;
};
