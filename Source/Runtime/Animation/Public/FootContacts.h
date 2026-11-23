#pragma once
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// 逐脚触地标签：每帧1字节，bit0=左脚、bit1=右脚。纯数据接口：不依赖存量头或GL，无共享可写状态。
//
// 原始判定（每脚独立、左右对称）：
//   (脚尖高 < ToeHeight 且 脚尖水平速度 < Speed) 或 (踝高 < AnkleHeight 且 踝3D速度 < Speed)
// 速度在片段内取中心差分，片段首尾取单侧差分，不跨片段接缝；单帧片段标签恒为0。
// 清理（片段内、顺序固定）：先把两侧均为触地、长度 <= MaxGap 的空隙填为触地，再删去长度 < MinContact 的触地段
// （含触及片段端点的段）。
inline constexpr uint8_t FootContactLeft = 1;
inline constexpr uint8_t FootContactRight = 2;

struct FFootContactThresholds
{
    float ToeHeight = 0.05f;
    float AnkleHeight = 0.11f;
    float Speed = 0.5f;
    int MaxGap = 2;
    int MinContact = 5;
    float StepSeconds = 0.0166667f;
};

struct FFootContactInput
{
    int Frames = 0;
    // 每帧依次为左踝、左趾、右踝、右趾的世界坐标xyz，共Frames*12个值
    std::vector<float> Joints;
    // 片段边界升序，首项0、末项Frames
    std::vector<int> Segments;
};

// 生成失败返回false并给出原因（输入维度、片段边界、阈值或坐标非有限）。
bool LabelFootContacts(const FFootContactInput& Input, const FFootContactThresholds& Thresholds,
                       std::vector<uint8_t>& OutLabels, std::string& Reason);

// FNV-1a 64位；Hash为续算的初值。
inline constexpr uint64_t Fnv1a64Offset = 14695981039346656037ull;
uint64_t HashBytesFnv1a64(const void* Data, size_t Size, uint64_t Hash = Fnv1a64Offset);
bool HashFileFnv1a64(const std::string& Path, uint64_t& OutHash, std::string& Reason);

// 文件格式（小端，无填充）：
//   0  char[8]  "MOOFC1\0\0"
//   8  u32      Version = 1
//   12 u32      Frames
//   16 i32[4]   骨索引：左踝、左趾、右踝、右趾
//   32 f32      ToeHeight, AnkleHeight, Speed
//   44 u32      MaxGap, MinContact
//   52 f32      StepSeconds
//   56 u64      源动作文件全字节FNV-1a64
//   64 u8[Frames] 标签（取值0..3），其后不得有尾随字节
inline constexpr uint32_t FootContactVersion = 1;
inline constexpr size_t FootContactHeaderSize = 64;
inline constexpr std::array<int32_t, 4> FootContactDefaultBones = {4, 5, 8, 9};

struct FFootContactFile
{
    std::array<int32_t, 4> Bones = FootContactDefaultBones;
    FFootContactThresholds Thresholds;
    uint64_t SourceHash = 0;
    std::vector<uint8_t> Labels;
};

enum class EFootContactError : uint8_t
{
    None, Open, Magic, Version, Truncated, FrameCount, Hash, Bones, Label, Trailing
};
const char* GetFootContactErrorName(EFootContactError Error);

std::vector<uint8_t> SerializeFootContacts(const FFootContactFile& File);
// ExpectedFrames/ExpectedHash为调用方由源动作得到的期望值；失败时Out不完整，Reason给出可读原因。
EFootContactError ParseFootContacts(const uint8_t* Data, size_t Size, int ExpectedFrames, uint64_t ExpectedHash,
                                    FFootContactFile& Out, std::string& Reason);
bool SaveFootContacts(const std::string& Path, const FFootContactFile& File, std::string& Reason);
EFootContactError LoadFootContacts(const std::string& Path, int ExpectedFrames, uint64_t ExpectedHash,
                                   FFootContactFile& Out, std::string& Reason);
