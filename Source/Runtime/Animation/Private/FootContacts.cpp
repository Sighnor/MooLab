#include "FootContacts.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iterator>

namespace
{
constexpr char FootContactMagic[8] = {'M', 'O', 'O', 'F', 'C', '1', '\0', '\0'};

// 显式按小端逐字节读写，不依赖主机字节序与结构体布局。
void PutU32(std::vector<uint8_t>& Out, uint32_t Value)
{
    for(int Byte = 0; Byte < 4; ++Byte) Out.push_back(static_cast<uint8_t>(Value >> (8 * Byte)));
}
void PutU64(std::vector<uint8_t>& Out, uint64_t Value)
{
    for(int Byte = 0; Byte < 8; ++Byte) Out.push_back(static_cast<uint8_t>(Value >> (8 * Byte)));
}
void PutF32(std::vector<uint8_t>& Out, float Value)
{
    uint32_t Bits = 0;
    std::memcpy(&Bits, &Value, sizeof(Bits));
    PutU32(Out, Bits);
}
uint32_t GetU32(const uint8_t* Data)
{
    uint32_t Value = 0;
    for(int Byte = 0; Byte < 4; ++Byte) Value |= static_cast<uint32_t>(Data[Byte]) << (8 * Byte);
    return Value;
}
uint64_t GetU64(const uint8_t* Data)
{
    uint64_t Value = 0;
    for(int Byte = 0; Byte < 8; ++Byte) Value |= static_cast<uint64_t>(Data[Byte]) << (8 * Byte);
    return Value;
}
float GetF32(const uint8_t* Data)
{
    const uint32_t Bits = GetU32(Data);
    float Value = 0.f;
    std::memcpy(&Value, &Bits, sizeof(Value));
    return Value;
}

// 片段内速度：中心差分，首尾单侧差分；Joint为0..3（左踝、左趾、右踝、右趾）。
void JointVelocity(const FFootContactInput& Input, int Frame, int Begin, int End, int Joint, float StepSeconds, double Out[3])
{
    const int Previous = Frame > Begin ? Frame - 1 : Frame;
    const int Next = Frame + 1 < End ? Frame + 1 : Frame;
    const double Seconds = (Next - Previous) * static_cast<double>(StepSeconds);
    const float* A = &Input.Joints[static_cast<size_t>(Previous) * 12 + Joint * 3];
    const float* B = &Input.Joints[static_cast<size_t>(Next) * 12 + Joint * 3];
    for(int Axis = 0; Axis < 3; ++Axis) Out[Axis] = (static_cast<double>(B[Axis]) - A[Axis]) / Seconds;
}

// 对[Begin, End)内单脚的0/1序列做两步清理。
void CleanRuns(std::vector<uint8_t>& Mask, int Begin, int End, int MaxGap, int MinContact)
{
    for(int Start = Begin; Start < End;)
    {
        int Stop = Start;
        while(Stop < End && Mask[Stop] == Mask[Start]) ++Stop;
        if(!Mask[Start] && Start > Begin && Stop < End && Stop - Start <= MaxGap)
            for(int Frame = Start; Frame < Stop; ++Frame) Mask[Frame] = 1;
        Start = Stop;
    }
    for(int Start = Begin; Start < End;)
    {
        int Stop = Start;
        while(Stop < End && Mask[Stop] == Mask[Start]) ++Stop;
        if(Mask[Start] && Stop - Start < MinContact)
            for(int Frame = Start; Frame < Stop; ++Frame) Mask[Frame] = 0;
        Start = Stop;
    }
}
}

bool LabelFootContacts(const FFootContactInput& Input, const FFootContactThresholds& Thresholds,
                       std::vector<uint8_t>& OutLabels, std::string& Reason)
{
    OutLabels.clear();
    if(Input.Frames <= 0 || Input.Joints.size() != static_cast<size_t>(Input.Frames) * 12)
    { Reason = "触地标注输入维度不符。"; return false; }
    const auto& Segments = Input.Segments;
    bool bSegments = Segments.size() >= 2 && Segments.front() == 0 && Segments.back() == Input.Frames;
    for(size_t Index = 1; bSegments && Index < Segments.size(); ++Index) bSegments = Segments[Index] > Segments[Index - 1];
    if(!bSegments) { Reason = "片段边界须严格升序且覆盖[0, Frames]。"; return false; }
    const bool bThresholds = std::isfinite(Thresholds.ToeHeight) && std::isfinite(Thresholds.AnkleHeight)
        && std::isfinite(Thresholds.Speed) && Thresholds.Speed > 0.f && std::isfinite(Thresholds.StepSeconds)
        && Thresholds.StepSeconds > 0.f && Thresholds.MaxGap >= 0 && Thresholds.MinContact >= 1;
    if(!bThresholds) { Reason = "触地阈值无效。"; return false; }
    for(float Value : Input.Joints)
        if(!std::isfinite(Value)) { Reason = "关节坐标含非有限值。"; return false; }

    OutLabels.assign(Input.Frames, 0);
    std::vector<uint8_t> Mask(Input.Frames, 0);
    for(int Side = 0; Side < 2; ++Side)
    {
        const int Ankle = Side * 2;
        const int Toe = Side * 2 + 1;
        std::fill(Mask.begin(), Mask.end(), static_cast<uint8_t>(0));
        for(size_t Segment = 0; Segment + 1 < Segments.size(); ++Segment)
        {
            const int Begin = Segments[Segment];
            const int End = Segments[Segment + 1];
            if(End - Begin < 2) continue;
            for(int Frame = Begin; Frame < End; ++Frame)
            {
                double ToeVelocity[3], AnkleVelocity[3];
                JointVelocity(Input, Frame, Begin, End, Toe, Thresholds.StepSeconds, ToeVelocity);
                JointVelocity(Input, Frame, Begin, End, Ankle, Thresholds.StepSeconds, AnkleVelocity);
                const float ToeHeight = Input.Joints[static_cast<size_t>(Frame) * 12 + Toe * 3 + 1];
                const float AnkleHeight = Input.Joints[static_cast<size_t>(Frame) * 12 + Ankle * 3 + 1];
                const double ToeSpeed = std::sqrt(ToeVelocity[0] * ToeVelocity[0] + ToeVelocity[2] * ToeVelocity[2]);
                const double AnkleSpeed = std::sqrt(AnkleVelocity[0] * AnkleVelocity[0] + AnkleVelocity[1] * AnkleVelocity[1]
                                                    + AnkleVelocity[2] * AnkleVelocity[2]);
                const bool bToe = ToeHeight < Thresholds.ToeHeight && ToeSpeed < Thresholds.Speed;
                const bool bAnkle = AnkleHeight < Thresholds.AnkleHeight && AnkleSpeed < Thresholds.Speed;
                Mask[Frame] = bToe || bAnkle ? 1 : 0;
            }
            CleanRuns(Mask, Begin, End, Thresholds.MaxGap, Thresholds.MinContact);
        }
        const uint8_t Bit = Side == 0 ? FootContactLeft : FootContactRight;
        for(int Frame = 0; Frame < Input.Frames; ++Frame) if(Mask[Frame]) OutLabels[Frame] |= Bit;
    }
    Reason.clear();
    return true;
}

uint64_t HashBytesFnv1a64(const void* Data, size_t Size, uint64_t Hash)
{
    const auto* Byte = static_cast<const unsigned char*>(Data);
    for(size_t Index = 0; Index < Size; ++Index) { Hash ^= Byte[Index]; Hash *= 1099511628211ull; }
    return Hash;
}

bool HashFileFnv1a64(const std::string& Path, uint64_t& OutHash, std::string& Reason)
{
    std::ifstream File(Path, std::ios::binary);
    if(!File) { Reason = "无法打开 " + Path; return false; }
    std::vector<char> Buffer(1 << 20);
    uint64_t Hash = Fnv1a64Offset;
    while(File)
    {
        File.read(Buffer.data(), static_cast<std::streamsize>(Buffer.size()));
        Hash = HashBytesFnv1a64(Buffer.data(), static_cast<size_t>(File.gcount()), Hash);
    }
    if(!File.eof()) { Reason = "读取失败 " + Path; return false; }
    OutHash = Hash;
    return true;
}

const char* GetFootContactErrorName(EFootContactError Error)
{
    switch(Error)
    {
    case EFootContactError::None: return "none";
    case EFootContactError::Open: return "open";
    case EFootContactError::Magic: return "magic";
    case EFootContactError::Version: return "version";
    case EFootContactError::Truncated: return "truncated";
    case EFootContactError::FrameCount: return "frame-count";
    case EFootContactError::Hash: return "hash";
    case EFootContactError::Bones: return "bones";
    case EFootContactError::Label: return "label";
    case EFootContactError::Trailing: return "trailing";
    }
    return "unknown";
}

std::vector<uint8_t> SerializeFootContacts(const FFootContactFile& File)
{
    std::vector<uint8_t> Out(std::begin(FootContactMagic), std::end(FootContactMagic));
    Out.reserve(FootContactHeaderSize + File.Labels.size());
    PutU32(Out, FootContactVersion);
    PutU32(Out, static_cast<uint32_t>(File.Labels.size()));
    for(int32_t Bone : File.Bones) PutU32(Out, static_cast<uint32_t>(Bone));
    PutF32(Out, File.Thresholds.ToeHeight);
    PutF32(Out, File.Thresholds.AnkleHeight);
    PutF32(Out, File.Thresholds.Speed);
    PutU32(Out, static_cast<uint32_t>(File.Thresholds.MaxGap));
    PutU32(Out, static_cast<uint32_t>(File.Thresholds.MinContact));
    PutF32(Out, File.Thresholds.StepSeconds);
    PutU64(Out, File.SourceHash);
    Out.insert(Out.end(), File.Labels.begin(), File.Labels.end());
    return Out;
}

EFootContactError ParseFootContacts(const uint8_t* Data, size_t Size, int ExpectedFrames, uint64_t ExpectedHash,
                                    FFootContactFile& Out, std::string& Reason)
{
    // 先头部、后长度、再内容：截断与尾随字节先于哈希/标签判定，保证错误类别唯一。
    if(Size < sizeof(FootContactMagic)) { Reason = "文件截断（不足magic）。"; return EFootContactError::Truncated; }
    if(std::memcmp(Data, FootContactMagic, sizeof(FootContactMagic)) != 0) { Reason = "magic不符。"; return EFootContactError::Magic; }
    if(Size < FootContactHeaderSize) { Reason = "文件截断（头部不完整）。"; return EFootContactError::Truncated; }
    const uint32_t Version = GetU32(Data + 8);
    if(Version != FootContactVersion) { Reason = "版本 " + std::to_string(Version) + " 不受支持。"; return EFootContactError::Version; }
    const uint32_t Frames = GetU32(Data + 12);
    if(ExpectedFrames < 0 || Frames != static_cast<uint32_t>(ExpectedFrames))
    {
        Reason = "帧数 " + std::to_string(Frames) + " 与源动作 " + std::to_string(ExpectedFrames) + " 不符。";
        return EFootContactError::FrameCount;
    }
    const size_t Expected = FootContactHeaderSize + Frames;
    if(Size < Expected) { Reason = "文件截断（标签不完整）。"; return EFootContactError::Truncated; }
    if(Size > Expected) { Reason = "文件含尾随字节。"; return EFootContactError::Trailing; }
    for(int Index = 0; Index < 4; ++Index) Out.Bones[Index] = static_cast<int32_t>(GetU32(Data + 16 + 4 * Index));
    if(Out.Bones != FootContactDefaultBones) { Reason = "骨索引与左踝4/左趾5/右踝8/右趾9约定不符。"; return EFootContactError::Bones; }
    Out.Thresholds.ToeHeight = GetF32(Data + 32);
    Out.Thresholds.AnkleHeight = GetF32(Data + 36);
    Out.Thresholds.Speed = GetF32(Data + 40);
    Out.Thresholds.MaxGap = static_cast<int>(GetU32(Data + 44));
    Out.Thresholds.MinContact = static_cast<int>(GetU32(Data + 48));
    Out.Thresholds.StepSeconds = GetF32(Data + 52);
    Out.SourceHash = GetU64(Data + 56);
    if(Out.SourceHash != ExpectedHash) { Reason = "源动作哈希不符：标签须由当前long_motion.bin重新生成。"; return EFootContactError::Hash; }
    Out.Labels.assign(Data + FootContactHeaderSize, Data + Expected);
    for(uint8_t Label : Out.Labels)
        if(Label > (FootContactLeft | FootContactRight)) { Reason = "标签值超出0..3。"; return EFootContactError::Label; }
    Reason.clear();
    return EFootContactError::None;
}

bool SaveFootContacts(const std::string& Path, const FFootContactFile& File, std::string& Reason)
{
    const std::vector<uint8_t> Bytes = SerializeFootContacts(File);
    std::ofstream Stream(Path, std::ios::binary | std::ios::trunc);
    if(!Stream) { Reason = "无法写入 " + Path; return false; }
    Stream.write(reinterpret_cast<const char*>(Bytes.data()), static_cast<std::streamsize>(Bytes.size()));
    Stream.close();
    if(!Stream) { Reason = "写入失败 " + Path; return false; }
    return true;
}

EFootContactError LoadFootContacts(const std::string& Path, int ExpectedFrames, uint64_t ExpectedHash,
                                   FFootContactFile& Out, std::string& Reason)
{
    std::ifstream Stream(Path, std::ios::binary);
    if(!Stream) { Reason = "无法打开 " + Path; return EFootContactError::Open; }
    const std::vector<uint8_t> Bytes((std::istreambuf_iterator<char>(Stream)), std::istreambuf_iterator<char>());
    if(Stream.bad()) { Reason = "读取失败 " + Path; return EFootContactError::Open; }
    return ParseFootContacts(Bytes.data(), Bytes.size(), ExpectedFrames, ExpectedHash, Out, Reason);
}
