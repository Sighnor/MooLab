// 仅由SandboxScene.cpp包含：Motion Matching私有实现。只读原始long_motion、database与逐脚触地标签，不写资源，
// 不调用会影响Physics基线的旧拼接/对齐函数。
namespace
{
constexpr int MatchBoneCount = 23;
constexpr int MatchRawFrames = 53500;
constexpr int MatchDatabaseRows = 53450;
constexpr int MatchLeftToe = 5;
constexpr int MatchRightToe = 9;
constexpr int MatchWindow = 75;
constexpr int MatchBlendSteps = 25;
constexpr int MatchMaxTracks = 5;
constexpr int MatchForcedSteps = 50;
constexpr int MatchDebounceSteps = 10;
constexpr int MatchColumnStride = 10;
constexpr int MatchSeams[] = {0, 281, 562, 13153, 25744, 39622, 53500};
constexpr float MatchPi = 3.14159265358979f;
constexpr float MatchMaxSpeed = 10.f;
constexpr float MatchMaxYawRate = 15.f;
constexpr float MatchBackwardSpeed = 0.5f;
constexpr float MatchRetargetSpeed = 0.5f;
constexpr float MatchVelocityTolerance = 0.1f;
constexpr float MatchFeetTolerance = 0.02f;
constexpr float MatchVelocityHalflife = 0.15f;
constexpr float MatchYawHalflife = 0.2f;

float WrapAngle(float Angle)
{
    return static_cast<float>(std::remainder(static_cast<double>(Angle), 2.0 * 3.14159265358979323846));
}
quat YawRotation(float Angle) { return quat(rad_to_deg(Angle), vec3(0.f, 1.f, 0.f)); }
vec3 Horizontal(const vec3& Value) { return vec3(Value.x, 0.f, Value.z); }
quat NormalizeRotation(const quat& Rotation)
{
    const float Norm = std::sqrt(dot(Rotation, Rotation));
    return Norm > 0.f ? Rotation * (1.f / Norm) : quat();
}
bool IsFinite(const vec3& Value) { return std::isfinite(Value.x) && std::isfinite(Value.y) && std::isfinite(Value.z); }
bool IsFinite(const quat& Rotation)
{
    return std::isfinite(Rotation.w) && std::isfinite(Rotation.x) && std::isfinite(Rotation.y) && std::isfinite(Rotation.z);
}
float SmoothStep(float Value) { return Value * Value * (3.f - 2.f * Value); }

// 播放轨：Frame为已显示的原始帧，End为窗口内最后可播放帧；根骨按独立偏航刚体变换对齐。
// Yaw为该轨根偏航的连续（不回绕）累积值，混合时按标量插值，避免轨间偏航差越过180度时最短弧翻转。
struct FMatchTrack
{
    int Frame = 0;
    int End = 0;
    int BlendStep = MatchBlendSteps;
    float Yaw = 0.f;
    quat Align;
    vec3 Offset;
};

class FMotionMatcher
{
public:
    // FootLabels为已校验的逐帧逐脚标签（bit0左、bit1右）；不足53500帧或取值越界时MM不可用。
    bool Initialize(const BVH_Motion& Source, const Database& Data, const std::vector<uint8_t>& FootLabels, std::string& Reason);
    bool IsAvailable() const { return Raw != nullptr && !Candidates.empty(); }
    int GetCandidateCount() const { return static_cast<int>(Candidates.size()); }
    void Reset();
    void Tick(const vec3& DesiredVelocity, const std::array<float, 3>& Weights);
    const array1d<vec3>& GetPositions() const { return Positions; }
    const array1d<quat>& GetRotations() const { return Rotations; }
    // 混合后的局部姿态；根骨为世界系，其余为父骨系偏移。
    const array1d<vec3>& GetLocalPositions() const { return LocalPositions; }
    // 当前显示帧的逐脚标签：取有效权重最大的轨（并列取较新轨）当前帧的标签。
    uint8_t GetFootContacts() const { return DisplayContacts; }
    int GetSearchQueryContacts() const { return SearchQueryContacts; }
    int GetSearchResultContacts() const { return SearchResultContacts; }
    bool IsSearchContactFallback() const { return bSearchContactFallback; }
    int GetSourceFrame() const { return TrackCount > 0 ? Tracks[TrackCount - 1].Frame : -1; }
    int GetClipStart() const { return ClipStart; }
    bool IsBlending() const { return TrackCount > 1; }
    int GetSearchStep() const { return SearchStep; }
    int GetSearchCount() const { return SearchCount; }
    float GetSearchCost() const { return SearchCost; }
    double GetSearchMs() const { return SearchMs; }
    int GetExcludedCount() const { return Excluded; }

private:
    int Search(const float* QueryVelocity, const float* QueryYawRate, const float* QueryFeet, int Contact,
               const std::array<float, 3>& Weights, float& BestCost) const;
    void StartClip(int Frame);
    void ComposeTrack(const FMatchTrack& Track, array1d<vec3>& OutPositions, array1d<quat>& OutRotations, quat& OutResidual) const;
    void ComposeBlend();
    void ForwardKinematics();
    void Advance();
    void RootFeatures(vec3& LeftToe, vec3& RightToe) const;
    float RootYaw() const { return Yaw(Rotations(0)); }

    const BVH_Motion* Raw = nullptr;
    std::vector<float> RawYaw;
    std::vector<int> Candidates;
    std::vector<uint8_t> FootContacts;
    std::vector<float> Velocities;   // 每帧6列水平局部速度(x,z)
    std::vector<float> YawRates;     // 每帧6列偏航角速
    std::vector<float> Feet;         // 每帧左右脚尖水平局部位置(x,z)
    int Excluded = 0;

    // 最旧轨为基底，较新轨依次以smoothstep权重淡入；新轨完全接管后丢弃其前各轨。
    std::array<FMatchTrack, MatchMaxTracks> Tracks;
    int TrackCount = 0;
    int ClipStart = -1;
    int StepIndex = 0;
    int StepsSinceSearch = 0;
    int SearchStep = -1;
    int SearchCount = 0;
    float SearchCost = 0.f;
    double SearchMs = 0.;
    int SearchQueryContacts = -1;
    int SearchResultContacts = -1;
    bool bSearchContactFallback = false;
    uint8_t DisplayContacts = 0;
    float TurnSign = 1.f;
    float BlendYaw = 0.f;
    vec3 SearchDesired;
    vec3 LastRoot;
    float LastYaw = 0.f;
    array1d<vec3> TrackPositions, LocalPositions, Positions;
    array1d<quat> TrackRotations, LocalRotations, Rotations;
};

bool FMotionMatcher::Initialize(const BVH_Motion& Source, const Database& Data, const std::vector<uint8_t>& FootLabels, std::string& Reason)
{
    Raw = nullptr;
    Candidates.clear();
    if(Source.nbones() != MatchBoneCount || Source.nframes() != MatchRawFrames || Source.bone_parents.size != MatchBoneCount)
    { Reason = "动作库维度不是53500x23。"; return false; }
    if(FootLabels.size() != static_cast<size_t>(MatchRawFrames)
        || std::any_of(FootLabels.begin(), FootLabels.end(), [](uint8_t Label) { return Label > 3; }))
    { Reason = "逐脚触地标签缺失、帧数不符或取值越界。"; return false; }
    bool bSkeleton = Source.bone_parents(0) == -1;
    for(int Bone = 1; Bone < MatchBoneCount; ++Bone)
        bSkeleton &= Source.bone_parents(Bone) >= 0 && Source.bone_parents(Bone) < Bone;
    bSkeleton &= Source.bone_parents(5) == 4 && Source.bone_parents(4) == 3 && Source.bone_parents(3) == 2
        && Source.bone_parents(2) == 1 && Source.bone_parents(9) == 8 && Source.bone_parents(8) == 7
        && Source.bone_parents(7) == 6 && Source.bone_parents(6) == 1;
    if(!bSkeleton) { Reason = "骨架父子索引与根骨0/LeftToe5/RightToe9约定不符。"; return false; }
    const int Rows = Data.nframes();
    if(Rows != MatchDatabaseRows || Data.contacts.size != Rows || Data.vels.rows != Rows || Data.vels.cols != 6
        || Data.avels.rows != Rows || Data.avels.cols != 6 || Data.features.rows != Rows || Data.features.cols < 2)
    { Reason = "database各数组维度不符。"; return false; }
    // database的单支撑脚标签只做维度/取值校验，不参与检索；检索使用逐脚标签。
    std::vector<int> RowOfFrame(Rows, -1);
    for(int Row = 0; Row < Rows; ++Row)
    {
        const int Frame = Data.frames(Row);
        if(Frame < 0 || Frame >= Rows || RowOfFrame[Frame] >= 0 || (Data.contacts(Row) != 0 && Data.contacts(Row) != 1))
        { Reason = "database帧映射或接触标签无效。"; return false; }
        RowOfFrame[Frame] = Row;
    }

    // 原始根偏航与逐帧最短角差；角速度特征在内存重算，不使用数据库avels。
    // 倒走帧：水平速度超过阈值且与根朝向夹角大于90度；含倒走帧的窗口不作候选，S键只能转身行走。
    RawYaw.assign(MatchRawFrames, 0.f);
    std::vector<double> CumulativeYaw(MatchRawFrames, 0.);
    std::vector<int> BadPrefix(MatchRawFrames + 1, 0);
    std::vector<int> BackwardPrefix(MatchRawFrames + 1, 0);
    for(int Frame = 0; Frame < MatchRawFrames; ++Frame)
    {
        RawYaw[Frame] = Yaw(Source.bone_local_rotations(Frame, 0));
        bool bBad = !std::isfinite(RawYaw[Frame]) || !IsFinite(Source.bone_local_positions(Frame, 0));
        bool bBackward = false;
        if(Frame > 0)
        {
            const float Delta = WrapAngle(RawYaw[Frame] - RawYaw[Frame - 1]);
            CumulativeYaw[Frame] = CumulativeYaw[Frame - 1] + Delta;
            const vec3 Step = Horizontal(Source.bone_local_positions(Frame, 0) - Source.bone_local_positions(Frame - 1, 0));
            bBad |= !(length(Step) / StepSeconds <= MatchMaxSpeed) || !(std::fabs(Delta) / StepSeconds <= MatchMaxYawRate);
            const vec3 Facing(std::sin(RawYaw[Frame]), 0.f, std::cos(RawYaw[Frame]));
            bBackward = length(Step) / StepSeconds > MatchBackwardSpeed && dot(Step, Facing) < 0.f;
        }
        BadPrefix[Frame + 1] = BadPrefix[Frame] + (bBad ? 1 : 0);
        BackwardPrefix[Frame + 1] = BackwardPrefix[Frame] + (bBackward ? 1 : 0);
    }

    Velocities.assign(static_cast<size_t>(MatchRawFrames) * 12, 0.f);
    YawRates.assign(static_cast<size_t>(MatchRawFrames) * 6, 0.f);
    Feet.assign(static_cast<size_t>(MatchRawFrames) * 4, 0.f);
    FootContacts = FootLabels;
    Excluded = 0;
    array1d<vec3> FkPositions(MatchBoneCount);
    array1d<quat> FkRotations(MatchBoneCount);
    for(int Segment = 0; Segment + 1 < static_cast<int>(std::size(MatchSeams)); ++Segment)
    {
        // 窗口[s-1, s+75]不得跨越片段接缝。
        for(int Frame = MatchSeams[Segment] + 1; Frame + MatchWindow < MatchSeams[Segment + 1]; ++Frame)
        {
            if(Frame >= Rows) break;
            const int Row = RowOfFrame[Frame];
            bool bValid = BadPrefix[Frame + MatchWindow + 1] - BadPrefix[Frame] == 0
                && BackwardPrefix[Frame + MatchWindow + 1] - BackwardPrefix[Frame] == 0;
            const quat Local = inv_quat(YawRotation(RawYaw[Frame]));
            float* Velocity = &Velocities[static_cast<size_t>(Frame) * 12];
            float* YawRate = &YawRates[static_cast<size_t>(Frame) * 6];
            for(int Column = 0; Column < 6; ++Column)
            {
                const vec3 Stored = Data.vels(Row, Column);
                const int Begin = Column == 0 ? Frame - 1 : Frame + MatchColumnStride * (Column - 1);
                const int End = Column == 0 ? Frame : Frame + MatchColumnStride * Column;
                const float Seconds = (End - Begin) * StepSeconds;
                const vec3 Recomputed = Local * Horizontal(Source.bone_local_positions(End, 0) - Source.bone_local_positions(Begin, 0)) / Seconds;
                Velocity[Column * 2 + 0] = Stored.x;
                Velocity[Column * 2 + 1] = Stored.z;
                YawRate[Column] = static_cast<float>((CumulativeYaw[End] - CumulativeYaw[Begin]) / Seconds);
                bValid &= IsFinite(Stored) && std::isfinite(YawRate[Column]) && length(Horizontal(Stored)) <= MatchMaxSpeed
                    && std::fabs(YawRate[Column]) <= MatchMaxYawRate
                    && length(Horizontal(Stored) - Recomputed) <= MatchVelocityTolerance;
            }
            batch_forward_kinematics_full(Source, Frame, FkPositions, FkRotations);
            float* Foot = &Feet[static_cast<size_t>(Frame) * 4];
            for(int Side = 0; Side < 2; ++Side)
            {
                const vec3 Stored = Data.features(Row, Side);
                const vec3 Recomputed = Local * (FkPositions(Side == 0 ? MatchLeftToe : MatchRightToe) - FkPositions(0));
                Foot[Side * 2 + 0] = Stored.x;
                Foot[Side * 2 + 1] = Stored.z;
                bValid &= IsFinite(Stored) && length(Horizontal(Stored - Recomputed)) <= MatchFeetTolerance;
            }
            if(bValid) Candidates.push_back(Frame);
            else ++Excluded;
        }
    }
    if(Candidates.empty()) { Reason = "没有满足窗口、接缝、倒走与有限性约束的候选帧。"; return false; }
    Raw = &Source;
    TrackPositions.resize(MatchBoneCount); LocalPositions.resize(MatchBoneCount); Positions.resize(MatchBoneCount);
    TrackRotations.resize(MatchBoneCount); LocalRotations.resize(MatchBoneCount); Rotations.resize(MatchBoneCount);
    Reason.clear();
    return true;
}

int FMotionMatcher::Search(const float* QueryVelocity, const float* QueryYawRate, const float* QueryFeet, int Contact,
                           const std::array<float, 3>& Weights, float& BestCost) const
{
    // 全量确定性argmin：严格小于且候选升序，并列时保留最小源帧。Contact>=0时只取逐脚标签完全相等的候选。
    int Best = -1;
    BestCost = std::numeric_limits<float>::infinity();
    for(int Frame : Candidates)
    {
        if(Contact >= 0 && FootContacts[Frame] != Contact) continue;
        const float* Velocity = &Velocities[static_cast<size_t>(Frame) * 12];
        const float* YawRate = &YawRates[static_cast<size_t>(Frame) * 6];
        const float* Foot = &Feet[static_cast<size_t>(Frame) * 4];
        float VelocityCost = 0.f, YawCost = 0.f, FeetCost = 0.f;
        for(int Index = 0; Index < 12; ++Index) VelocityCost += (QueryVelocity[Index] - Velocity[Index]) * (QueryVelocity[Index] - Velocity[Index]);
        for(int Index = 0; Index < 6; ++Index) YawCost += (QueryYawRate[Index] - YawRate[Index]) * (QueryYawRate[Index] - YawRate[Index]);
        for(int Index = 0; Index < 4; ++Index) FeetCost += (QueryFeet[Index] - Foot[Index]) * (QueryFeet[Index] - Foot[Index]);
        const float Cost = Weights[0] * VelocityCost / 6.f + Weights[1] * YawCost / (6.f * 1.6f * 1.6f)
            + Weights[2] * FeetCost / (2.f * 0.25f * 0.25f);
        if(std::isfinite(Cost) && Cost < BestCost) { BestCost = Cost; Best = Frame; }
    }
    return Best;
}

void FMotionMatcher::ComposeTrack(const FMatchTrack& Track, array1d<vec3>& OutPositions, array1d<quat>& OutRotations, quat& OutResidual) const
{
    for(int Bone = 1; Bone < MatchBoneCount; ++Bone)
    {
        OutPositions(Bone) = Raw->bone_local_positions(Track.Frame, Bone);
        OutRotations(Bone) = Raw->bone_local_rotations(Track.Frame, Bone);
    }
    const vec3 Root = Raw->bone_local_positions(Track.Frame, 0);
    OutPositions(0) = Track.Offset + Track.Align * Horizontal(Root) + vec3(0.f, Root.y, 0.f);
    // 根旋转拆为偏航与非偏航残差；残差接近单位，可安全slerp。
    const quat Rotation = NormalizeRotation(Track.Align * Raw->bone_local_rotations(Track.Frame, 0));
    OutResidual = NormalizeRotation(inv_quat(YawRotation(Yaw(Rotation))) * Rotation);
}

void FMotionMatcher::ComposeBlend()
{
    quat Residual, TrackResidual;
    ComposeTrack(Tracks[0], LocalPositions, LocalRotations, Residual);
    BlendYaw = Tracks[0].Yaw;
    std::array<float, MatchMaxTracks> Weights{};
    Weights[0] = 1.f;
    for(int Index = 1; Index < TrackCount; ++Index)
    {
        ComposeTrack(Tracks[Index], TrackPositions, TrackRotations, TrackResidual);
        const float Weight = SmoothStep(static_cast<float>(Tracks[Index].BlendStep) / MatchBlendSteps);
        Weights[Index] = Weight;
        BlendYaw += (Tracks[Index].Yaw - BlendYaw) * Weight;
        Residual = NormalizeRotation(slerp(Residual, TrackResidual, Weight));
        LocalPositions(0) = lerp(LocalPositions(0), TrackPositions(0), Weight);
        for(int Bone = 1; Bone < MatchBoneCount; ++Bone)
        {
            LocalPositions(Bone) = lerp(LocalPositions(Bone), TrackPositions(Bone), Weight);
            LocalRotations(Bone) = NormalizeRotation(slerp(LocalRotations(Bone), TrackRotations(Bone), Weight));
        }
    }
    LocalRotations(0) = NormalizeRotation(YawRotation(WrapAngle(BlendYaw)) * Residual);
    // 第i轨有效权重 w_i * Π_{j>i}(1 - w_j)；由新到旧遍历并取严格更大者，并列时保留较新轨。
    float Remaining = 1.f, BestWeight = -1.f;
    for(int Index = TrackCount - 1; Index >= 0; --Index)
    {
        const float Effective = Weights[Index] * Remaining;
        if(Effective > BestWeight) { BestWeight = Effective; DisplayContacts = FootContacts[Tracks[Index].Frame]; }
        Remaining *= 1.f - Weights[Index];
    }
}

void FMotionMatcher::ForwardKinematics()
{
    for(int Bone = 0; Bone < MatchBoneCount; ++Bone)
    {
        const int Parent = Raw->bone_parents(Bone);
        if(Parent < 0) { Rotations(Bone) = LocalRotations(Bone); Positions(Bone) = LocalPositions(Bone); continue; }
        Rotations(Bone) = Rotations(Parent) * LocalRotations(Bone);
        Positions(Bone) = Positions(Parent) + Rotations(Parent) * LocalPositions(Bone);
    }
}

void FMotionMatcher::StartClip(int Frame)
{
    // 查询描述当前显示帧，匹配帧s对齐到当前显示根位置/朝向，下一步播放s+1。
    // 新轨权重从0淡入、旧轨照常播放，因此混合中途再次搜索也不会使姿态跳变。
    if(TrackCount == MatchMaxTracks)
    {
        std::move(Tracks.begin() + 1, Tracks.end(), Tracks.begin());
        --TrackCount;
    }
    FMatchTrack& Track = Tracks[TrackCount];
    Track.Frame = Frame;
    Track.End = Frame + MatchWindow;
    Track.BlendStep = TrackCount == 0 ? MatchBlendSteps : 0;
    Track.Yaw = BlendYaw;
    Track.Align = YawRotation(WrapAngle(BlendYaw - RawYaw[Frame]));
    Track.Offset = Horizontal(Positions(0)) - Track.Align * Horizontal(Raw->bone_local_positions(Frame, 0));
    ++TrackCount;
    ClipStart = Frame;
}

void FMotionMatcher::Advance()
{
    for(int Index = 0; Index < TrackCount; ++Index)
    {
        FMatchTrack& Track = Tracks[Index];
        const int Next = std::min(Track.Frame + 1, Track.End);
        Track.Yaw += WrapAngle(RawYaw[Next] - RawYaw[Track.Frame]);
        Track.Frame = Next;
        if(Index > 0) Track.BlendStep = std::min(Track.BlendStep + 1, MatchBlendSteps);
    }
    // 各轨同移2π整数倍，保持轨间差不变并限制连续偏航的量级。
    const float Shift = Tracks[0].Yaw - WrapAngle(Tracks[0].Yaw);
    for(int Index = 0; Index < TrackCount; ++Index) Tracks[Index].Yaw -= Shift;
    ComposeBlend();
    for(int Index = TrackCount - 1; Index > 0; --Index)
    {
        if(Tracks[Index].BlendStep < MatchBlendSteps) continue;
        std::move(Tracks.begin() + Index, Tracks.begin() + TrackCount, Tracks.begin());
        TrackCount -= Index;
        break;
    }
    ForwardKinematics();
}

void FMotionMatcher::RootFeatures(vec3& LeftToe, vec3& RightToe) const
{
    const quat Local = inv_quat(YawRotation(RootYaw()));
    LeftToe = Local * (Positions(MatchLeftToe) - Positions(0));
    RightToe = Local * (Positions(MatchRightToe) - Positions(0));
}

void FMotionMatcher::Reset()
{
    MOO_CHECK_ALWAYS(IsAvailable(), "Motion Matching reset without valid data");
    StepIndex = 0;
    StepsSinceSearch = 0;
    SearchCount = 0;
    SearchMs = 0.;
    TurnSign = 1.f;
    SearchDesired = vec3(0.f);
    ClipStart = -1;
    TrackCount = 0;
    BlendYaw = 0.f;
    // 首次搜索静止片段：零速度/零角速，不使用脚部项与接触过滤。
    const float Zero[12] = {};
    const auto Begin = std::chrono::steady_clock::now();
    const int Frame = Search(Zero, Zero, Zero, -1, {1.f, 1.f, 0.f}, SearchCost);
    MOO_CHECK_ALWAYS(Frame >= 0, "Motion Matching initial search found no candidate");
    SearchMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - Begin).count();
    SearchStep = 0;
    SearchCount = 1;
    SearchQueryContacts = -1;
    SearchResultContacts = FootContacts[Frame];
    bSearchContactFallback = false;
    Positions.zero();
    Rotations.set(quat());
    StartClip(Frame);
    ComposeBlend();
    ForwardKinematics();
    LastRoot = Positions(0);
    LastYaw = RootYaw();
}

void FMotionMatcher::Tick(const vec3& DesiredVelocity, const std::array<float, 3>& Weights)
{
    MOO_CHECK_ALWAYS(Raw != nullptr && TrackCount > 0, "Motion Matching ticked before Reset");
    const vec3 Root = Positions(0);
    const float CurrentYaw = RootYaw();
    const quat Local = inv_quat(YawRotation(CurrentYaw));
    const vec3 Velocity = Local * Horizontal(Root - LastRoot) / StepSeconds;
    const float YawRate = WrapAngle(CurrentYaw - LastYaw) / StepSeconds;
    const vec3 Desired = Local * Horizontal(DesiredVelocity);
    const float DesiredSpeed = length(Desired);
    const float MeasuredSpeed = length(Velocity);
    // 局部系朝向为+z。朝向转角=atan2(cross,dot)最短角；近180度沿用上次转向，初始为正。
    const float BaseHeading = MeasuredSpeed > 0.1f ? std::atan2(Velocity.x, Velocity.z) : 0.f;
    float FacingTurn = 0.f, VelocityTurn = 0.f;
    if(DesiredSpeed > 1e-4f)
    {
        const float TargetHeading = std::atan2(Desired.x, Desired.z);
        FacingTurn = TargetHeading;
        // 近180度带内延续上次转向：符号相反时取绕远路的等价转角，而非镜像取值。
        if(std::fabs(FacingTurn) > MatchPi - 0.1f) { if(FacingTurn * TurnSign < 0.f) FacingTurn += TurnSign * 2.f * MatchPi; }
        else if(std::fabs(FacingTurn) > 1e-3f) TurnSign = FacingTurn > 0.f ? 1.f : -1.f;
        VelocityTurn = WrapAngle(TargetHeading - BaseHeading);
        if(std::fabs(VelocityTurn) > MatchPi - 0.1f && VelocityTurn * TurnSign < 0.f) VelocityTurn += TurnSign * 2.f * MatchPi;
    }
    // 预测：速率按一阶半衰期趋向期望；速度方向与朝向按临界阻尼弹簧转过完整转角，
    // 反向指令沿弧线转身而非经零速倒走。零输入时期望速度/转角为零，角速衰减到零。
    const float VelocityRate = std::log(2.f) / MatchVelocityHalflife;
    const float YawDamping = 2.f * std::log(2.f) / MatchYawHalflife;
    const float YawSlope = YawDamping * FacingTurn - YawRate;
    auto YawError = [&](float Time) { return (FacingTurn + YawSlope * Time) * std::exp(-YawDamping * Time); };
    auto HeadingRemaining = [&](float Time) { return (1.f + YawDamping * Time) * std::exp(-YawDamping * Time); };
    vec3 LeftToe, RightToe;
    RootFeatures(LeftToe, RightToe);
    float QueryVelocity[12], QueryYawRate[6], QueryFeet[4] = {LeftToe.x, LeftToe.z, RightToe.x, RightToe.z};
    QueryVelocity[0] = Velocity.x;
    QueryVelocity[1] = Velocity.z;
    QueryYawRate[0] = YawRate;
    for(int Column = 1; Column < 6; ++Column)
    {
        // 第k列对应区间[10(k-1),10k]步：速度取区间中点tau=(k-0.5)*10*dt，角速取区间精确平均。
        const float Begin = (Column - 1) * MatchColumnStride * StepSeconds;
        const float End = Column * MatchColumnStride * StepSeconds;
        const float Middle = 0.5f * (Begin + End);
        const float Speed = DesiredSpeed + (MeasuredSpeed - DesiredSpeed) * std::exp(-VelocityRate * Middle);
        const float Heading = BaseHeading + VelocityTurn * (1.f - HeadingRemaining(Middle));
        QueryVelocity[Column * 2 + 0] = Speed * std::sin(Heading);
        QueryVelocity[Column * 2 + 1] = Speed * std::cos(Heading);
        QueryYawRate[Column] = (YawError(Begin) - YawError(End)) / (End - Begin);
    }
    // 查询状态为当前显示帧的逐脚标签（比较前已由ComposeBlend选定）。
    const int Contact = DisplayContacts;

    ++StepsSinceSearch;
    FMatchTrack& Current = Tracks[TrackCount - 1];
    const vec3 PreviousDesired = Horizontal(SearchDesired);
    const bool bMoving = length(Horizontal(DesiredVelocity)) > 0.1f;
    const bool bWasMoving = length(PreviousDesired) > 0.1f;
    bool bChanged = bMoving != bWasMoving || length(Horizontal(DesiredVelocity) - PreviousDesired) > 0.5f;
    if(bMoving && bWasMoving)
        bChanged |= dot(normalize(Horizontal(DesiredVelocity)), normalize(PreviousDesired)) < std::cos(deg_to_rad(30.f));
    // 实际速度与期望速度持续偏离超过0.5m/s时同样按防抖间隔重搜，避免停步/转身时沿用不合适片段直到强制重搜。
    bChanged |= length(Velocity - Desired) > MatchRetargetSpeed;
    const bool bSearch = StepsSinceSearch >= MatchForcedSteps || Current.Frame >= Current.End
        || (StepsSinceSearch >= MatchDebounceSteps && bChanged);
    bool bFinite = true;
    for(float Value : QueryVelocity) bFinite &= std::isfinite(Value);
    for(float Value : QueryYawRate) bFinite &= std::isfinite(Value);
    for(float Value : QueryFeet) bFinite &= std::isfinite(Value);
    if(bSearch && bFinite)
    {
        const auto Begin = std::chrono::steady_clock::now();
        float Cost = 0.f;
        int Best = Search(QueryVelocity, QueryYawRate, QueryFeet, Contact, Weights, Cost);
        bSearchContactFallback = Best < 0;
        if(Best < 0) Best = Search(QueryVelocity, QueryYawRate, QueryFeet, -1, Weights, Cost);
        SearchMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - Begin).count();
        SearchStep = StepIndex;
        ++SearchCount;
        StepsSinceSearch = 0;
        SearchDesired = DesiredVelocity;
        SearchQueryContacts = Contact;
        SearchResultContacts = Best >= 0 ? FootContacts[Best] : -1;
        if(Best >= 0)
        {
            SearchCost = Cost;
            // 最优即当前轨的当前/下一帧时只延长窗口，不做自身过渡。
            if(Best == Current.Frame || Best == Current.Frame + 1) { Current.End = Best + MatchWindow; ClipStart = Best; }
            else StartClip(Best);
        }
    }
    LastRoot = Root;
    LastYaw = CurrentYaw;
    Advance();
    ++StepIndex;
}
}
