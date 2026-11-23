// 仅由SandboxScene.cpp包含（位于MotionMatching.inl之后）：MM显示姿态的脚锁IK后处理。
// 只改两腿UpLeg/Leg/Foot全局旋转与膝、踝、趾位置；趾保持动画全局旋转，不回馈Matcher，不触碰Physics。
//
// 锁定：标签上升沿记录当前显示脚尖的水平位置L，y沿用动画；绕脚尖向上抬跟实现soft伸展距离。
// 膝角优先且每步最多14.5度；抬跟仍不可达时如实报告残差，不直接将脚尖位置写成锁点。
// 释放条件仍为关闭/标签下降/水平距离>0.3m/残差>2cm；世界系脚尖偏移与soft权重共同smoothstep淡出。
// N至少8步，按偏移或抬跟弧长/0.01m确定；N结束后仍受膝角带约束，直到膝/踝/趾与动画重合后交还。
// 重合判据与显示姿态同为float世界坐标、以米计；权重归零后至多FootIkMaxSettleSteps步必定交还。
// 新的上升沿可打断淡出并在当前显示脚尖处重新锁定。
namespace
{
constexpr int FootIkChains[2][4] = {{2, 3, 4, 5}, {6, 7, 8, 9}};
constexpr float FootIkReleaseDistance = 0.3f;
constexpr double FootIkReleaseClamp = 0.02;
constexpr int FootIkMinFadeSteps = 8;
constexpr float FootIkFadeStepLength = 0.01f;
constexpr int FootIkCalibrationStride = 7;
constexpr float FootIkCalibrationMinSine = 0.2f;
// 权重归零后膝角带每步14.5度，13步即可覆盖180度；30步仍未重合视为判据失效，强制交还并计数（探针要求为0）。
constexpr int FootIkMaxSettleSteps = 30;

// 交还容差（m）。双精度重解的输入是float世界坐标，误差约为坐标ULP（≈FLT_EPSILON*|x|）
// 经近伸直膝的几何放大（约10倍）；32*FLT_EPSILON*|x|覆盖该放大，1e-5m下限覆盖原点附近，且远小于可见位移。
float FootIkHandbackTolerance(const vec3& Position)
{
    const float Scale = std::max({1.f, std::fabs(Position.x), std::fabs(Position.y), std::fabs(Position.z)});
    return std::max(1e-5f, 32.f * std::numeric_limits<float>::epsilon() * Scale);
}

std::array<double, 3> ToDouble(const vec3& Value) { return {Value.x, Value.y, Value.z}; }
std::array<double, 4> ToDouble(const quat& Value) { return {Value.w, Value.x, Value.y, Value.z}; }
vec3 ToVec(const std::array<double, 3>& Value)
{
    return vec3(static_cast<float>(Value[0]), static_cast<float>(Value[1]), static_cast<float>(Value[2]));
}
quat ToQuat(const std::array<double, 4>& Value)
{
    return quat(static_cast<float>(Value[0]), static_cast<float>(Value[1]), static_cast<float>(Value[2]), static_cast<float>(Value[3]));
}

struct FFootIkLeg
{
    bool bLabel = false;
    bool bLocked = false;
    bool bFading = false;
    vec3 Lock;
    int FadeStep = 0;
    int FadeSteps = 0;
    vec3 FadeDelta;
    double PreviousKnee = -1.;
    vec3 LiftAxis = vec3(1.f, 0.f, 0.f);
    float LiftAngle = 0.f;
    bool bConflict = false;
    int ConflictSteps = 0;
    int FadeConflictSteps = 0;
    int InvalidSteps = 0;
    // 权重归零后已等待的步数、观察到的最大等待步数、超限强制交还次数。
    int SettleSteps = 0;
    int MaxSettleSteps = 0;
    int ForcedHandbacks = 0;
    float ClampError = 0.f;
    vec3 AnimToe;
    EFootRelease LastRelease = EFootRelease::None;
    std::array<int, FootReleaseKinds> ReleaseCounts{};
};

class FFootLocker
{
public:
    // 由原始动作标定两膝在UpLeg局部系的弯曲参考轴；只读源动作。
    void Initialize(const BVH_Motion& Source);
    void Reset() { Legs = {}; LastMs = 0.; }
    // Contacts为当前显示帧的逐脚标签；Local为Matcher混合后的局部偏移；Positions/Rotations为待后处理的全局姿态。
    void Apply(uint8_t Contacts, bool bEnabled, const array1d<vec3>& Local, array1d<vec3>& Positions, array1d<quat>& Rotations);
    const FFootIkLeg& GetLeg(int Side) const { return Legs[Side]; }
    // 标定样本与平均轴的平均点积（1为完全一致），供验证报告。
    float GetHingeConsistency(int Side) const { return HingeConsistency[Side]; }
    double GetLastMs() const { return LastMs; }

private:
    struct FLegPose
    {
        vec3 Knee, Ankle, Toe;
        quat Upper, Lower, Foot;
        vec3 LiftAxis;
        double Lift = 0.;
        double KneeAngle = 0.;
        double Unreachable = 0.;
        bool bConflict = false;
        bool bSolved = false;
    };
    FLegPose SolveLeg(int Side, const vec3& ToeTarget, double SoftWeight, const array1d<vec3>& Local,
                      const array1d<vec3>& Positions, const array1d<quat>& Rotations) const;

    std::array<vec3, 2> HingeAxes = {vec3(1.f, 0.f, 0.f), vec3(1.f, 0.f, 0.f)};
    std::array<float, 2> HingeConsistency = {0.f, 0.f};
    std::array<FFootIkLeg, 2> Legs;
    double LastMs = 0.;
};

void FFootLocker::Initialize(const BVH_Motion& Source)
{
    array1d<vec3> FkPositions(Source.nbones());
    array1d<quat> FkRotations(Source.nbones());
    std::array<std::vector<vec3>, 2> Samples;
    for(int Frame = 0; Frame < Source.nframes(); Frame += FootIkCalibrationStride)
    {
        batch_forward_kinematics_full(Source, Frame, FkPositions, FkRotations);
        for(int Side = 0; Side < 2; ++Side)
        {
            const vec3 Hip = FkPositions(FootIkChains[Side][0]);
            const vec3 Reach = FkPositions(FootIkChains[Side][2]) - Hip;
            const vec3 Thigh = FkPositions(FootIkChains[Side][1]) - Hip;
            const vec3 Normal = cross(Reach, Thigh);
            const float Scale = length(Reach) * length(Thigh);
            if(!(Scale > 0.f) || length(Normal) <= FootIkCalibrationMinSine * Scale) continue;
            Samples[Side].push_back(inv_quat(FkRotations(FootIkChains[Side][0])) * (Normal / length(Normal)));
        }
    }
    for(int Side = 0; Side < 2; ++Side)
    {
        vec3 Sum(0.f);
        for(const vec3& Sample : Samples[Side]) Sum = Sum + Sample;
        if(Samples[Side].empty() || !(length(Sum) > 0.f)) continue;
        HingeAxes[Side] = Sum / length(Sum);
        double Agreement = 0.;
        for(const vec3& Sample : Samples[Side]) Agreement += dot(Sample, HingeAxes[Side]);
        HingeConsistency[Side] = static_cast<float>(Agreement / Samples[Side].size());
    }
}

FFootLocker::FLegPose FFootLocker::SolveLeg(int Side, const vec3& ToeTarget, double SoftWeight,
                                            const array1d<vec3>& Local, const array1d<vec3>& Positions, const array1d<quat>& Rotations) const
{
    const int* Chain = FootIkChains[Side];
    FFootPivotIkInput Input;
    Input.Chain.Root = ToDouble(Positions(Chain[0]));
    Input.Chain.UpperRotation = ToDouble(Rotations(Chain[0]));
    Input.Chain.LowerRotation = ToDouble(Rotations(Chain[1]));
    Input.Chain.MidOffset = ToDouble(Local(Chain[1]));
    Input.Chain.EndOffset = ToDouble(Local(Chain[2]));
    Input.Chain.ReferenceAxis = ToDouble(Rotations(Chain[0]) * HingeAxes[Side]);
    Input.FootRotation = ToDouble(Rotations(Chain[2]));
    Input.ToeOffset = ToDouble(Local(Chain[3]));
    Input.ToeTarget = ToDouble(ToeTarget);
    Input.PreviousLiftAxis = ToDouble(Legs[Side].LiftAxis);
    Input.PreviousKneeAngle = Legs[Side].PreviousKnee;
    Input.SoftWeight = SoftWeight;
    FFootPivotIkOutput Output;
    FLegPose Pose;
    Pose.bSolved = SolveFootPivotIk(Input, Output);
    if(!Pose.bSolved)
    {
        // 非有限输入：保持动画姿态
        Pose.Knee = Positions(Chain[1]);
        Pose.Ankle = Positions(Chain[2]);
        Pose.Toe = Positions(Chain[3]);
        Pose.Upper = Rotations(Chain[0]);
        Pose.Lower = Rotations(Chain[1]);
        Pose.Foot = Rotations(Chain[2]);
        return Pose;
    }
    Pose.Upper = ToQuat(Output.Chain.UpperRotation);
    Pose.Lower = ToQuat(Output.Chain.LowerRotation);
    Pose.Foot = ToQuat(Output.FootRotation);
    Pose.Knee = ToVec(Output.Chain.Mid);
    Pose.Ankle = ToVec(Output.Chain.End);
    Pose.Toe = ToVec(Output.Toe);
    Pose.LiftAxis = ToVec(Output.LiftAxis);
    Pose.Lift = Output.LiftAngle;
    Pose.KneeAngle = Output.KneeAngle;
    Pose.Unreachable = Output.Residual;
    Pose.bConflict = Output.bConflict;
    return Pose;
}

void FFootLocker::Apply(uint8_t Contacts, bool bEnabled, const array1d<vec3>& Local, array1d<vec3>& Positions, array1d<quat>& Rotations)
{
    const auto Begin = std::chrono::steady_clock::now();
    for(int Side = 0; Side < 2; ++Side)
    {
        FFootIkLeg& Leg = Legs[Side];
        const int* Chain = FootIkChains[Side];
        const bool bLabel = (Contacts & (Side == 0 ? 1 : 2)) != 0;
        const bool bRising = bLabel && !Leg.bLabel;
        const vec3 AnimAnkle = Positions(Chain[2]);
        const vec3 AnimToe = Positions(Chain[3]);
        const vec3 KneeToHip = Positions(Chain[0]) - Positions(Chain[1]);
        const vec3 KneeToAnkle = AnimAnkle - Positions(Chain[1]);
        const double AnimKnee = std::atan2(static_cast<double>(length(cross(KneeToHip, KneeToAnkle))),
                                          static_cast<double>(dot(KneeToHip, KneeToAnkle)));
        auto LockTarget = [&]() { return vec3(Leg.Lock.x, AnimToe.y, Leg.Lock.z); };
        auto Solve = [&](const vec3& Target, double Weight)
        {
            FLegPose Result = SolveLeg(Side, Target, Weight, Local, Positions, Rotations);
            Leg.InvalidSteps += !Result.bSolved; // 每次调用都计数，包括淡出被新上升沿打断的情况。
            return Result;
        };
        FLegPose Pose;
        bool bApply = false;
        EFootRelease Release = EFootRelease::None;
        if(Leg.bLocked)
        {
            if(!bEnabled) Release = EFootRelease::Disabled;
            else if(!bLabel) Release = EFootRelease::Label;
            else if(length(Horizontal(Leg.Lock - AnimToe)) > FootIkReleaseDistance) Release = EFootRelease::Distance;
            Pose = Solve(LockTarget(), 1.);
            bApply = true;
            if(Release == EFootRelease::None && Pose.Unreachable > FootIkReleaseClamp) Release = EFootRelease::Clamp;
        }
        else
        {
            if(Leg.bFading)
            {
                Leg.FadeStep = std::min(Leg.FadeStep + 1, Leg.FadeSteps);
                const float Weight = 1.f - SmoothStep(static_cast<float>(Leg.FadeStep) / Leg.FadeSteps);
                // 淡出保留完整抬跟可行区间；将旧抬跟角作为下限会在动画转脚时强迫不可达目标。
                Pose = Solve(AnimToe + Leg.FadeDelta * Weight, Weight);
                bApply = true;
                if(Weight == 0.f && Pose.bSolved)
                {
                    ++Leg.SettleSteps;
                    Leg.MaxSettleSteps = std::max(Leg.MaxSettleSteps, Leg.SettleSteps);
                    const float Tolerance = FootIkHandbackTolerance(AnimToe);
                    const bool bMatched = !Pose.bConflict && length(Pose.Knee - Positions(Chain[1])) <= Tolerance
                        && length(Pose.Ankle - AnimAnkle) <= Tolerance && length(Pose.Toe - AnimToe) <= Tolerance;
                    if(bMatched || Leg.SettleSteps >= FootIkMaxSettleSteps)
                    {
                        Leg.ForcedHandbacks += !bMatched;
                        Leg.SettleSteps = 0;
                        Leg.bFading = false;
                        bApply = false; // 交还原始动画，不留数值重解残差。
                    }
                }
            }
            if(bRising && bEnabled)
            {
                Leg.Lock = bApply ? Pose.Toe : AnimToe;
                Leg.bLocked = true;
                Leg.bFading = false;
                Leg.SettleSteps = 0;
                Pose = Solve(LockTarget(), 1.);
                bApply = true;
                if(Pose.Unreachable > FootIkReleaseClamp) Release = EFootRelease::Clamp;
            }
        }
        if(bApply && !Pose.bSolved)
        {
            // 无效输入无法安全延续IK：回到本步动画并清状态；探针要求InvalidSteps为0，不隐瞒回退。
            Leg.bLocked = Leg.bFading = false;
            Leg.SettleSteps = 0;
            bApply = false;
            Release = EFootRelease::None;
        }
        Leg.ClampError = bApply ? static_cast<float>(Pose.Unreachable) : 0.f;
        Leg.LiftAngle = bApply ? static_cast<float>(Pose.Lift) : 0.f;
        Leg.bConflict = bApply && Pose.bConflict;
        Leg.ConflictSteps += Leg.bConflict;
        Leg.FadeConflictSteps += Leg.bConflict && Leg.bFading;
        if(Release != EFootRelease::None)
        {
            Leg.bLocked = false;
            Leg.bFading = true;
            Leg.FadeStep = 0;
            Leg.SettleSteps = 0;
            Leg.FadeDelta = Pose.Toe - AnimToe;
            const double FadeTravel = std::max(static_cast<double>(length(Leg.FadeDelta)), length(Local(Chain[3])) * Pose.Lift);
            Leg.FadeSteps = std::max(FootIkMinFadeSteps, static_cast<int>(std::ceil(FadeTravel / FootIkFadeStepLength)));
            Leg.LastRelease = Release;
            ++Leg.ReleaseCounts[static_cast<int>(Release)];
        }
        if(bApply)
        {
            Rotations(Chain[0]) = Pose.Upper;
            Rotations(Chain[1]) = Pose.Lower;
            Rotations(Chain[2]) = Pose.Foot;
            Positions(Chain[1]) = Pose.Knee;
            Positions(Chain[2]) = Pose.Ankle;
            Positions(Chain[3]) = Pose.Toe;
        }
        Leg.PreviousKnee = bApply ? Pose.KneeAngle : AnimKnee;
        if(bApply) Leg.LiftAxis = Pose.LiftAxis;
        Leg.bLabel = bLabel;
        Leg.AnimToe = AnimToe;
    }
    LastMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - Begin).count();
}
}
