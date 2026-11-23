#include "SandboxScene.h"
#include "GpuRenderer.h"
#include "PbrShading.h"
#include "FootContacts.h"
#include "TwoBoneIk.h"
#include "Core.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <sstream>
#include <tuple>
#include <type_traits>
#include <vector>

// 唯一消费旧头定义的编译单元。私有inl不参与独立编译。
#include "character.hpp"
#include "controller.hpp"
#include "database.hpp"
#include "motion.hpp"
#include "curve.hpp"
#include "renderer.hpp"
#include "simulator.hpp"

namespace
{
constexpr float StepSeconds = 0.0166667f;
// CpuRenderer.models 与 GpuItems 的模型序号：光源、角色，其后为按骨骼顺序的胶囊。
constexpr int LightModelIndex = 0;
constexpr int CharacterModelIndex = 1;
constexpr int FirstCapsuleModelIndex = 2;
void CopyRenderMatrix(float Out[16], const mat4& Matrix)
{
    static_assert(sizeof(mat4) == 16 * sizeof(float) && std::is_trivially_copyable_v<mat4>);
    std::memcpy(Out, &Matrix, sizeof(mat4));
}
FMeshDataView MeshView(const MooMesh& Mesh)
{
    return {Mesh.vertices, Mesh.normals, Mesh.texcoords, Mesh.indices, Mesh.vertex_count, Mesh.triangle_count};
}
FImageRGBFView ImageView(const texture& Image)
{
    static_assert(sizeof(vec3) == 3 * sizeof(float));
    return {reinterpret_cast<const float*>(Image.colors), Image.cols, Image.rows};
}
std::array<float, 3> VecView(const vec3& Value) { return {Value.x, Value.y, Value.z}; }
float Yaw(const quat& Rotation)
{
    const vec3 Forward = Rotation * vec3(0.f, 0.f, 1.f);
    return std::atan2(Forward.x, Forward.z);
}
bool InRange(double Value, double Min, double Max)
{
    return std::isfinite(Value) && Value >= Min && Value <= Max;
}
bool IsValidMaterial(const FSurfaceMaterial& Material)
{
    return InRange(Material.Metallic, 0., 1.) && InRange(Material.Roughness, 0.05, 1.);
}
// 资源按编译期绝对目录定位，不切换进程工作目录；目录缺失时警告并退回工作目录相对的旧路径。
std::string ResourceRoot()
{
#ifdef MOOLAB_RESOURCE_DIR
    std::error_code Error;
    if(std::filesystem::is_directory(MOOLAB_RESOURCE_DIR, Error)) return std::string(MOOLAB_RESOURCE_DIR) + "/";
    std::cerr << "Warning: resource dir not found, falling back to ../resources: " << MOOLAB_RESOURCE_DIR << '\n';
#endif
    return "../resources/";
}
struct FGpuDrawItem
{
    MooModel* Model = nullptr;
    int ModelIndex = 0;
    FMeshId Mesh;
    EShading Shading = EShading::Light;
    FPbrParams Pbr;
};
}

#include "MotionMatching.inl"
#include "FootIk.inl"

namespace
{
constexpr int FootContactSegments[] = {0, 281, 562, 13153, 25744, 39622, 53500};

// 读取并校验逐脚标签：帧数与源动作一致、哈希为long_motion.bin全字节FNV-1a64。
bool LoadMatchingFootContacts(const std::string& Path, const std::string& MotionPath, int Frames, std::vector<uint8_t>& Labels,
                              std::string& Reason)
{
    uint64_t Hash = 0;
    if(!HashFileFnv1a64(MotionPath, Hash, Reason)) return false;
    FFootContactFile File;
    const EFootContactError Error = LoadFootContacts(Path, Frames, Hash, File, Reason);
    if(Error != EFootContactError::None)
    {
        Reason = std::string(GetFootContactErrorName(Error)) + "：" + Reason;
        return false;
    }
    Labels = std::move(File.Labels);
    return true;
}
}

struct FSandboxScene::FImpl
{
    FSceneSettings Settings;
    ESceneMode Mode = ESceneMode::Physics;
    bool bPaused = false;
    bool bStep = false;
    int NextFrame = 0;
    int LastFrame = -1;
    std::string Message;
    MooCamera Camera{1, 0.1f, 0.1f, 0.05f, 50.f, vec3(0.f, 1.3f, -6.f), mat3(), 720, 720};
    Camera_Controller CameraController{};
    Character CharacterData;
    Database DatabaseData;
    BVH_Motion RawMotion, Motion, ActiveMotion;
    Simulator Physics, MatchShapes;
    FMotionMatcher Matcher;
    FFootLocker FootLocker;
    std::string MatchingReason;
    FSceneInput LastValidInput;
    Point_Light Light{vec3(0.f, 1.f, 1.5f), vec3(0.f, 0.f, -1.f), vec3(20.f)};
    MooMesh SkinMesh;
    array1d<MooMesh> CapsuleMeshes;
    MooModel LightModel, SkinModel;
    array1d<MooModel> CapsuleModels;
    updated_paramters Parameters{};
    MooRenderer CpuRenderer;
    array1d<vec3> Positions, LastPositions;
    array1d<quat> Rotations, LastRotations;
    std::vector<FGpuDrawItem> GpuItems;

    FImpl(const std::string& FootContactsPath);
    void ResetView();
    void ResetPhysics();
    void ResetMatching();
    void ApplyFootIk();
    void TickPhysics();
    void TickMatching(const FSceneInput& Input);
    const Simulator& Shapes() const { return Mode == ESceneMode::MotionMatching ? MatchShapes : Physics; }
    void UpdatePresentation(const FSceneInput& Input);
    bool ValidateSettings(const FSceneSettings& Value);
    // 把Settings中逐对象的材质同步到CPU模型与GPU绘制项；只改标量，不触碰纹理所有权。
    void ApplyMaterials();
    const FSurfaceMaterial* FindModelMaterial(int ModelIndex) const;
    // 当前模式实际使用的地面平面（原始系数）；每次按模式计算，不缓存。
    std::array<double, 4> GetGroundPlane() const;
    FFrameParams MakeFrameParams();
    bool MakeGroundParams(const FFrameParams& Frame, int Width, int Height, FGroundParams& Out) const;
    // CPU地面pass：清屏后、模型前逐像素写入地面颜色与view depth。
    void DrawGroundCpu(const FFrameParams& FrameParams, FBO& Target);
    void DrawContacts(FGpuRenderer* Renderer);
};

FSandboxScene::FImpl::FImpl(const std::string& FootContactsPath)
{
    const std::string Resources = ResourceRoot();
    auto Resource = [&Resources](const char* Name) { return Resources + Name; };
    Camera.set_view(vec3(0.f, 0.f, 1.f));
    Camera.update_orientation(vec3(0.f, 1.f, 0.f));
    bind_controller(CameraController, &Camera.position, &Camera.direction);
    character_load(CharacterData, Resource("character.bin").c_str());
    Database_load(DatabaseData, Resource("database.bin").c_str());
    Motion_load(RawMotion, Resource("long_motion.bin").c_str());
    // 保持历史拼接的顺序、常数和浮点运算；只移除未消费的实验缓冲。原始动作另供Motion Matching只读使用。
    BVH_Motion Motion1, Motion2, Motion3, Motion4, Motion5, Motion6;
    Motion1 = motion_sub_sequence(RawMotion, 0, 281);
    Motion2 = motion_sub_sequence(RawMotion, 281, 562);
    Motion3 = motion_sub_sequence(RawMotion, 562, 13153);
    Motion4 = motion_sub_sequence(RawMotion, 13153, 25744);
    Motion5 = motion_sub_sequence(RawMotion, 25744, 39622);
    Motion6 = motion_sub_sequence(RawMotion, 39622, 53500);
    Motion = concatenate_two_motions(Motion1, Motion2, 281, 30);
    Motion = concatenate_two_motions(Motion, Motion3, 562, 30);
    Motion = concatenate_two_motions(Motion, Motion4, 13153, 30);
    Motion = concatenate_two_motions(Motion, Motion5, 25744, 30);
    Motion = concatenate_two_motions(Motion, Motion6, 39622, 30);
    bind_simulator(Physics, Motion, 12600, 0.05);
    bind_simulator(MatchShapes, Motion, 12600, 0.05);
    MOO_CHECK_ALWAYS(Physics.bone_shapes.size == static_cast<int>(std::tuple_size_v<decltype(FSceneSnapshot::Bones)>)
                     && MatchShapes.bone_shapes.size == Physics.bone_shapes.size,
                     "capsule count %d does not match scene snapshot", Physics.bone_shapes.size);
    SkinMesh = make_character_rest_mesh(CharacterData);
    CapsuleMeshes.resize(Physics.bone_masses.size);
    CapsuleModels.resize(Physics.bone_masses.size);
    for(int Index = 0; Index < CapsuleMeshes.size; ++Index)
    {
        CapsuleMeshes(Index) = capsule(Physics.bone_shapes(Index).radius, Physics.bone_shapes(Index).length,
                                      18, 18, Physics.bone_shapes(Index).radius);
    }
    // 材质标量由ApplyMaterials按Settings逐对象写入；albedo纹理为各模型浅拷贝共享。
    MooMaterial Material(PBR);
    Material.tex = img_to_tex(Resource("skybox.png").c_str());
    LightModel = get_light_model(Light);
    SkinModel = mesh_material_to_model(SkinMesh, Material);
    for(int Index = 0; Index < CapsuleModels.size; ++Index)
        CapsuleModels(Index) = mesh_material_to_model(CapsuleMeshes(Index), Material);
    Parameters.light_dir = &Light.light_direction;
    Parameters.light_pos = &Light.light_position;
    Parameters.light_radiance = &Light.light_radiance;
    Parameters.shadow_map = Light.shadow_map;
    CpuRenderer.models.push_back(&LightModel);
    CpuRenderer.models.push_back(&SkinModel);
    for(int Index = 0; Index < CapsuleModels.size; ++Index) CpuRenderer.models.push_back(&CapsuleModels(Index));
    MOO_CHECK_ALWAYS(CpuRenderer.models[LightModelIndex] == &LightModel && CpuRenderer.models[CharacterModelIndex] == &SkinModel
                     && CpuRenderer.models.size() == static_cast<size_t>(FirstCapsuleModelIndex + SceneCapsuleCount),
                     "unexpected scene model order");
    ApplyMaterials();
    CpuRenderer.point_lights.push_back(&Light);
    Positions.resize(Motion.nbones());
    Rotations.resize(Motion.nbones());
    LastPositions.resize(Motion.nbones());
    LastRotations.resize(Motion.nbones());
    ActiveMotion = motion_sub_sequence(Motion, 12600, 13600);
    std::vector<uint8_t> FootLabels;
    const auto LabelBegin = std::chrono::steady_clock::now();
    const std::string LabelPath = FootContactsPath.empty() ? Resource("foot_contacts.bin") : FootContactsPath;
    std::string LabelReason;
    if(!LoadMatchingFootContacts(LabelPath, Resource("long_motion.bin"), RawMotion.nframes(), FootLabels, LabelReason))
        MatchingReason = "逐脚触地标签无效（" + LabelReason + "），请用FootContactBake重新生成。";
    else
    {
        std::cout << "Foot contacts load+hash: "
                  << std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - LabelBegin).count() << " ms\n";
        if(Matcher.Initialize(RawMotion, DatabaseData, FootLabels, MatchingReason)) FootLocker.Initialize(RawMotion);
    }
    if(!Matcher.IsAvailable()) std::cerr << "Motion Matching unavailable: " << MatchingReason << '\n';
    ResetPhysics();
}

void FSandboxScene::FImpl::ResetView()
{
    Camera.position = vec3(0.f, 1.3f, -6.f);
    Camera.set_view(vec3(0.f, 0.f, 1.f));
    Camera.update_orientation(vec3(0.f, 1.f, 0.f));
    CameraController = Camera_Controller{};
    bind_controller(CameraController, &Camera.position, &Camera.direction);
    Light.light_position = vec3(0.f, 1.f, 1.5f);
    Light.light_radiance = vec3(20.f);
    NextFrame = 0;
    LastFrame = -1;
    bStep = false;
    LastValidInput = FSceneInput{};
}

void FSandboxScene::FImpl::ResetPhysics()
{
    // bind只重写动态数组；不重建模型、材质、FBO或GPU登记指针。
    bind_simulator(Physics, Motion, 12600, 0.05);
    batch_forward_kinematics_full(Motion, 0, Positions, Rotations);
    LastPositions = Positions;
    LastRotations = Rotations;
    ResetView();
}

void FSandboxScene::FImpl::ResetMatching()
{
    // 独立的MM姿态与胶囊映射；不触碰Physics求解状态。首帧已触地的脚即刻锁定。
    Matcher.Reset();
    FootLocker.Reset();
    Positions = Matcher.GetPositions();
    Rotations = Matcher.GetRotations();
    ApplyFootIk();
    LastPositions = Positions;
    LastRotations = Rotations;
    update_simulator(MatchShapes, LastPositions, LastRotations, Positions, Rotations);
    ResetView();
    UpdatePresentation(FSceneInput{});
}

void FSandboxScene::FImpl::TickMatching(const FSceneInput& Input)
{
    FSceneInput Move = Input;
    if(!std::isfinite(Move.MoveX) || !std::isfinite(Move.MoveY)) Move = LastValidInput;
    else LastValidInput = Move;
    // WASD/左摇杆按相机水平朝向换算为世界期望速度，模长饱和到速度上限。
    // 相机近竖直时水平分量退化，退回世界+z而不归一化零向量。
    const vec3 CameraForward = Horizontal(Camera.direction);
    const float CameraForwardLength = length(CameraForward);
    const vec3 Forward = CameraForwardLength > 1e-4f ? CameraForward / CameraForwardLength : vec3(0.f, 0.f, 1.f);
    const vec3 Right = cross(Forward, vec3(0.f, 1.f, 0.f));
    vec3 Desired = Forward * (-Move.MoveY) + Right * Move.MoveX;
    const float Magnitude = length(Desired);
    if(Magnitude > 1.f) Desired = Desired / Magnitude;
    Matcher.Tick(Desired * Settings.MoveSpeed, Settings.MatchingWeights);
    Positions = Matcher.GetPositions();
    Rotations = Matcher.GetRotations();
    ApplyFootIk();
    update_simulator(MatchShapes, LastPositions, LastRotations, Positions, Rotations);
}

void FSandboxScene::FImpl::ApplyFootIk()
{
    // 只后处理显示姿态；Matcher内部姿态与下一步查询不受影响。
    FootLocker.Apply(Matcher.GetFootContacts(), Settings.bFootIk, Matcher.GetLocalPositions(), Positions, Rotations);
}

void FSandboxScene::FImpl::TickPhysics()
{
    for(int Substep = 0; Substep < Settings.Substeps; ++Substep)
    {
        Physics.simulate_gravity(Settings.Gravity);
        Physics.simulate_damp(Settings.Damping);
        const auto& Plane = Settings.ContactPlane;
        Physics.simulate_contact(Plane[0], Plane[1], Plane[2], Plane[3], Settings.ContactEpsilon);
        if(NextFrame < Settings.ControlEndFrame)
        {
            const int SourceFrame = std::min(static_cast<int>(0.5f * NextFrame), ActiveMotion.nframes() - 1);
            if(Settings.Control == EControlMode::Global)
                Physics.simulate_global_control(ActiveMotion.bone_local_rotations(SourceFrame));
            else if(Settings.Control == EControlMode::Local)
                Physics.simulate_local_control(ActiveMotion.bone_local_rotations(SourceFrame));
        }
        Physics.simulate(StepSeconds / Settings.Substeps);
    }
    deform_character_anim_bones(Physics, Positions, Rotations);
}

void FSandboxScene::FImpl::UpdatePresentation(const FSceneInput& Input)
{
    deform_character_anim_mesh(CharacterData, Positions, Rotations, SkinMesh);
    const Simulator& Source = Shapes();
    for(int Index = 0; Index < CapsuleModels.size; ++Index)
        CapsuleModels(Index).transform = mat4(eye3(), Source.bone_shapes(Index).pos)
            * mat4(quat_to_Rodrigues(Source.bone_shapes(Index).rot), vec3());
    // 非有限的视角输入不进入相机控制器，避免角度状态永久变为NaN。
    const bool bFiniteLook = std::isfinite(Input.LookX) && std::isfinite(Input.LookY);
    CameraController.dir_gamepad_control(bFiniteLook ? vec3(Input.LookX, 0.f, Input.LookY) : vec3(0.f), Settings.CameraSensitivity);
    Camera.update_orientation(vec3(0.f, 1.f, 0.f));
    CameraController.pos_gamepad_control(vec3(Positions(Bone_Hips).x, Settings.CameraHeight, Positions(Bone_Hips).z),
                                        Camera.orientation, Settings.CameraDistance);
    Light.light_position = Settings.bFollowLight ? Positions(Bone_Hips) + Camera.orientation * vec3(0.f, 0.f, 5.f)
        : vec3(Settings.LightPosition[0], Settings.LightPosition[1], Settings.LightPosition[2]);
    Light.light_radiance = Settings.LightBase + Settings.LightAmplitude * sin(deg_to_rad(0.5f * NextFrame));
    LightModel.transform = mat4(eye3(), Light.light_position) * mat4(0.01f);
}

bool FSandboxScene::FImpl::ValidateSettings(const FSceneSettings& Value)
{
    const auto& Plane = Value.ContactPlane;
    bool bValid = InRange(Value.Gravity, 0., 30.) && InRange(Value.Damping, 0., 5.)
        && Value.Substeps >= 1 && Value.Substeps <= 8 && Value.ControlEndFrame >= 0
        && Value.ControlEndFrame <= 1000000 && Value.Control <= EControlMode::Off
        && InRange(Plane[0], -1., 1.) && InRange(Plane[1], -1., 1.) && InRange(Plane[2], -1., 1.)
        && Plane[0] * Plane[0] + Plane[1] * Plane[1] + Plane[2] * Plane[2] > 0.
        && InRange(Plane[3], -2., 2.) && InRange(Value.ContactEpsilon, 0., 0.04)
        && InRange(Value.CameraDistance, 0.2, 20.) && InRange(Value.CameraHeight, -5., 5.)
        && InRange(Value.CameraSensitivity, 0., 0.2) && InRange(Value.LightBase, 0., 100.)
        && InRange(Value.LightAmplitude, 0., Value.LightBase) && InRange(Value.MoveSpeed, 0., 10.)
        && IsValidMaterial(Value.CharacterMaterial) && IsValidMaterial(Value.Ground.Material)
        && InRange(Value.Ground.CellSize, 0.05, 10.) && InRange(Value.Ground.MatchingHeight, -2., 2.);
    for(const FSurfaceMaterial& Material : Value.CapsuleMaterials) bValid &= IsValidMaterial(Material);
    double WeightSum = 0.;
    for(float Weight : Value.MatchingWeights) { bValid &= InRange(Weight, 0., 10.); WeightSum += Weight; }
    for(float Coordinate : Value.LightPosition) bValid &= InRange(Coordinate, -100., 100.);
    for(int Weight : Value.LegacyWeights) bValid &= Weight >= 0 && Weight <= 50;
    bValid &= WeightSum > 0.;
    if(!bValid) Message = "参数被拒绝：数值须有限、在允许范围内，平面法向非零且匹配权重不能全零。";
    else Message.clear();
    return bValid;
}

const FSurfaceMaterial* FSandboxScene::FImpl::FindModelMaterial(int ModelIndex) const
{
    if(ModelIndex == CharacterModelIndex) return &Settings.CharacterMaterial;
    const int Capsule = ModelIndex - FirstCapsuleModelIndex;
    if(Capsule >= 0 && Capsule < SceneCapsuleCount) return &Settings.CapsuleMaterials[Capsule];
    return nullptr;
}

void FSandboxScene::FImpl::ApplyMaterials()
{
    for(size_t Index = 0; Index < CpuRenderer.models.size(); ++Index)
    {
        const FSurfaceMaterial* Material = FindModelMaterial(static_cast<int>(Index));
        if(!Material) continue;
        MooModel& Model = *CpuRenderer.models[Index];
        for(int Slot = 0; Slot < Model.material_count; ++Slot)
        {
            Model.materials[Slot].metallic = Material->Metallic;
            Model.materials[Slot].roughness = Material->Roughness;
        }
    }
    for(auto& Item : GpuItems)
    {
        const FSurfaceMaterial* Material = FindModelMaterial(Item.ModelIndex);
        if(!Material || Item.Shading != EShading::Pbr) continue;
        Item.Pbr.Metallic = Material->Metallic;
        Item.Pbr.Roughness = Material->Roughness;
    }
}

std::array<double, 4> FSandboxScene::FImpl::GetGroundPlane() const
{
    if(Mode == ESceneMode::MotionMatching) return {0., 1., 0., -static_cast<double>(Settings.Ground.MatchingHeight)};
    return Settings.ContactPlane;
}

FSandboxScene::FSandboxScene(const std::string& FootContactsPath) : Impl(std::make_unique<FImpl>(FootContactsPath)) {}
FSandboxScene::~FSandboxScene() = default;
const FSceneSettings& FSandboxScene::GetSettings() const { return Impl->Settings; }
bool FSandboxScene::SetSettings(const FSceneSettings& Settings)
{
    if(!Impl->ValidateSettings(Settings)) return false;
    Impl->Settings = Settings;
    Impl->ApplyMaterials();
    return true;
}
bool FSandboxScene::ApplyCommand(const FSceneCommand& Command)
{
    if(Command.bSetMode && Command.Mode != ESceneMode::Physics && Command.Mode != ESceneMode::MotionMatching)
    {
        Impl->Message = "模式被拒绝：未知场景模式。";
        return false;
    }
    if(Command.bSetMode && Command.Mode == ESceneMode::MotionMatching && !Impl->Matcher.IsAvailable())
    {
        Impl->Message = "Motion Matching 不可用：" + Impl->MatchingReason;
        return false;
    }
    if(Command.bRestoreDefaults)
    {
        const bool bDefaultsValid = SetSettings(FSceneSettings{});
        MOO_CHECK_ALWAYS(bDefaultsValid, "default scene settings failed validation");
    }
    if(Command.bSetMode) Impl->Mode = Command.Mode;
    // 切模式本身即初始化目标模式；与重置同帧下达时只重置一次。
    if(Command.bSetMode || Command.bReset) Reset();
    if(Command.bTogglePause) Impl->bPaused = !Impl->bPaused;
    if(Command.bSingleStep) { Impl->bPaused = true; Impl->bStep = true; }
    return true;
}
void FSandboxScene::Reset()
{
    if(Impl->Mode == ESceneMode::MotionMatching) Impl->ResetMatching();
    else Impl->ResetPhysics();
}
void FSandboxScene::Tick(const FSceneInput& Input)
{
    if(Impl->bPaused && !Impl->bStep) return;
    Impl->bStep = false;
    Impl->LastPositions = Impl->Positions;
    Impl->LastRotations = Impl->Rotations;
    if(Impl->Mode == ESceneMode::MotionMatching) Impl->TickMatching(Input);
    else Impl->TickPhysics();
    Impl->UpdatePresentation(Input);
    Impl->LastFrame = Impl->NextFrame++;
}
FSceneSnapshot FSandboxScene::GetSnapshot() const
{
    FSceneSnapshot Result;
    Result.Mode = Impl->Mode;
    Result.bPaused = Impl->bPaused;
    Result.NextFrame = Impl->NextFrame;
    Result.LastFrame = Impl->LastFrame;
    Result.Message = Impl->Message;
    Result.RootPosition = VecView(Impl->Positions(0));
    Result.RootYaw = Yaw(Impl->Rotations(0));
    const vec3 Delta = Impl->Positions(0) - Impl->LastPositions(0);
    Result.RootSpeed = std::sqrt(Delta.x * Delta.x + Delta.z * Delta.z) / StepSeconds;
    Result.CameraPosition = VecView(Impl->Camera.position);
    Result.CameraYaw = Impl->CameraController.ang.x;
    Result.LightPosition = VecView(Impl->Light.light_position);
    Result.LightRadiance = VecView(Impl->Light.light_radiance);
    for(int Index = 0; Index < static_cast<int>(Result.Bones.size()); ++Index)
    {
        const auto& Shape = Impl->Shapes().bone_shapes(Index);
        Result.Bones[Index].Position = VecView(Shape.pos);
        Result.Bones[Index].Rotation = {Shape.rot.w, Shape.rot.x, Shape.rot.y, Shape.rot.z};
        Result.Bones[Index].bContact = Shape.ifcontact;
    }
    Result.GroundPlane = Impl->GetGroundPlane();
    const auto& Matcher = Impl->Matcher;
    Result.bMatchingAvailable = Matcher.IsAvailable();
    Result.CandidateCount = Matcher.GetCandidateCount();
    Result.ExcludedCount = Matcher.GetExcludedCount();
    if(Impl->Mode == ESceneMode::MotionMatching)
    {
        Result.SourceFrame = Matcher.GetSourceFrame();
        Result.ClipStart = Matcher.GetClipStart();
        Result.bBlending = Matcher.IsBlending();
        Result.SearchFrame = Matcher.GetSearchStep();
        Result.SearchCount = Matcher.GetSearchCount();
        Result.SearchCost = Matcher.GetSearchCost();
        Result.SearchMs = Matcher.GetSearchMs();
        Result.FootContacts = Matcher.GetFootContacts();
        Result.SearchQueryContacts = Matcher.GetSearchQueryContacts();
        Result.SearchResultContacts = Matcher.GetSearchResultContacts();
        Result.bSearchContactFallback = Matcher.IsSearchContactFallback();
        Result.FootIkMs = Impl->FootLocker.GetLastMs();
        for(int Side = 0; Side < 2; ++Side)
        {
            const FFootIkLeg& Leg = Impl->FootLocker.GetLeg(Side);
            FFootView& View = Result.Feet[Side];
            View.HipPosition = VecView(Impl->Positions(FootIkChains[Side][0]));
            View.KneePosition = VecView(Impl->Positions(FootIkChains[Side][1]));
            View.AnklePosition = VecView(Impl->Positions(FootIkChains[Side][2]));
            View.ToePosition = VecView(Impl->Positions(FootIkChains[Side][3]));
            View.AnimToePosition = VecView(Leg.AnimToe);
            View.LockPosition = VecView(Leg.Lock);
            View.bContactLabel = Leg.bLabel;
            View.bLocked = Leg.bLocked;
            View.bFading = Leg.bFading;
            View.ClampError = Leg.ClampError;
            View.LiftAngle = Leg.LiftAngle;
            View.bConflict = Leg.bConflict;
            View.ConflictSteps = Leg.ConflictSteps;
            View.FadeConflictSteps = Leg.FadeConflictSteps;
            View.InvalidSteps = Leg.InvalidSteps;
            View.MaxSettleSteps = Leg.MaxSettleSteps;
            View.ForcedHandbacks = Leg.ForcedHandbacks;
            View.LastRelease = Leg.LastRelease;
            View.ReleaseCounts = Leg.ReleaseCounts;
        }
    }
    return Result;
}

std::string GetSandboxResourcePath(const char* Name) { return ResourceRoot() + Name; }

bool BuildFootContactsFile(std::vector<uint8_t>& OutBytes, FFootContactStats& Stats, std::string& Report)
{
    OutBytes.clear();
    Stats = FFootContactStats{};
    const std::string MotionPath = GetSandboxResourcePath("long_motion.bin");
    FFootContactFile File;
    if(!HashFileFnv1a64(MotionPath, File.SourceHash, Report)) return false;
    BVH_Motion Source;
    Motion_load(Source, MotionPath.c_str());
    const int Frames = Source.nframes();
    if(Frames != FootContactSegments[std::size(FootContactSegments) - 1] || Source.nbones() <= 9)
    { Report = "long_motion.bin维度与片段边界不符。"; return false; }
    bool bSkeleton = true;
    for(int Side = 0; Side < 2; ++Side)
        for(int Link = 1; Link < 4; ++Link)
            bSkeleton &= Source.bone_parents(FootIkChains[Side][Link]) == FootIkChains[Side][Link - 1];
    if(!bSkeleton) { Report = "骨架与左右腿链2-3-4-5/6-7-8-9约定不符。"; return false; }

    FFootContactInput Input;
    Input.Frames = Frames;
    Input.Joints.resize(static_cast<size_t>(Frames) * 12);
    Input.Segments.assign(std::begin(FootContactSegments), std::end(FootContactSegments));
    std::vector<float> Root(static_cast<size_t>(Frames) * 2);
    array1d<vec3> FkPositions(Source.nbones());
    array1d<quat> FkRotations(Source.nbones());
    for(int Frame = 0; Frame < Frames; ++Frame)
    {
        batch_forward_kinematics_full(Source, Frame, FkPositions, FkRotations);
        for(int Joint = 0; Joint < 4; ++Joint)
        {
            const vec3 Value = FkPositions(File.Bones[Joint]);
            float* Out = &Input.Joints[static_cast<size_t>(Frame) * 12 + Joint * 3];
            Out[0] = Value.x;
            Out[1] = Value.y;
            Out[2] = Value.z;
        }
        Root[Frame * 2 + 0] = FkPositions(0).x;
        Root[Frame * 2 + 1] = FkPositions(0).z;
    }
    if(!LabelFootContacts(Input, File.Thresholds, File.Labels, Report)) return false;
    OutBytes = SerializeFootContacts(File);

    // 统计（只读）
    const auto& Labels = File.Labels;
    Stats.Frames = Frames;
    int Left = 0, Right = 0, Both = 0, Flight = 0;
    for(uint8_t Label : Labels)
    {
        Left += (Label & FootContactLeft) != 0;
        Right += (Label & FootContactRight) != 0;
        Both += Label == 3;
        Flight += Label == 0;
    }
    Stats.ContactRate = {static_cast<double>(Left) / Frames, static_cast<double>(Right) / Frames};
    Stats.BothRate = static_cast<double>(Both) / Frames;
    Stats.FlightRate = static_cast<double>(Flight) / Frames;
    for(int Pair = 0; Pair < 3; ++Pair)
    {
        const int A = FootContactSegments[Pair * 2], B = FootContactSegments[Pair * 2 + 1], End = FootContactSegments[Pair * 2 + 2];
        if(B - A != End - B) continue;
        int Agree = 0;
        for(int Offset = 0; Offset < B - A; ++Offset)
        {
            const uint8_t First = Labels[A + Offset], Second = Labels[B + Offset];
            Agree += ((First & 1) != 0) == ((Second & 2) != 0) && ((First & 2) != 0) == ((Second & 1) != 0);
        }
        Stats.MirrorAgreement[Pair] = static_cast<double>(Agree) / (B - A);
    }
    int StillBoth = 0;
    for(size_t Segment = 0; Segment + 1 < std::size(FootContactSegments); ++Segment)
    {
        const int Begin = FootContactSegments[Segment], End = FootContactSegments[Segment + 1];
        for(int Frame = Begin; Frame < End; ++Frame)
        {
            const int Previous = Frame > Begin ? Frame - 1 : Frame, Next = Frame + 1 < End ? Frame + 1 : Frame;
            if(Next == Previous) continue;
            const double Seconds = (Next - Previous) * static_cast<double>(File.Thresholds.StepSeconds);
            const double VX = (static_cast<double>(Root[Next * 2]) - Root[Previous * 2]) / Seconds;
            const double VZ = (static_cast<double>(Root[Next * 2 + 1]) - Root[Previous * 2 + 1]) / Seconds;
            if(std::sqrt(VX * VX + VZ * VZ) >= 0.1) continue;
            ++Stats.StillFrames;
            StillBoth += Labels[Frame] == 3;
        }
    }
    Stats.StillBothRate = Stats.StillFrames > 0 ? static_cast<double>(StillBoth) / Stats.StillFrames : 0.;
    Database Data;
    Database_load(Data, GetSandboxResourcePath("database.bin").c_str());
    int Supported = 0, Rows = 0;
    for(int Row = 0; Row < Data.nframes(); ++Row)
    {
        const int Frame = Data.frames(Row);
        if(Frame < 0 || Frame >= Frames) continue;
        ++Rows;
        Supported += (Labels[Frame] & (Data.contacts(Row) == 0 ? FootContactLeft : FootContactRight)) != 0;
    }
    Stats.SupportAgreement = Rows > 0 ? static_cast<double>(Supported) / Rows : 0.;
    std::vector<int> Runs;
    for(size_t Segment = 0; Segment + 1 < std::size(FootContactSegments); ++Segment)
        for(uint8_t Bit : {FootContactLeft, FootContactRight})
            for(int Frame = FootContactSegments[Segment]; Frame < FootContactSegments[Segment + 1];)
            {
                int Stop = Frame;
                while(Stop < FootContactSegments[Segment + 1] && ((Labels[Stop] & Bit) != 0) == ((Labels[Frame] & Bit) != 0)) ++Stop;
                if(Labels[Frame] & Bit) Runs.push_back(Stop - Frame);
                Frame = Stop;
            }
    std::sort(Runs.begin(), Runs.end());
    Stats.MinRun = Runs.empty() ? 0 : Runs.front();
    Stats.MedianRun = Runs.empty() ? 0 : Runs[Runs.size() / 2];

    std::ostringstream Text;
    Text << "frames=" << Frames << " sourceHash=" << std::hex << File.SourceHash << std::dec
         << " contactRate L=" << Stats.ContactRate[0] << " R=" << Stats.ContactRate[1] << " both=" << Stats.BothRate
         << " flight=" << Stats.FlightRate << " mirror=" << Stats.MirrorAgreement[0] << '/' << Stats.MirrorAgreement[1] << '/'
         << Stats.MirrorAgreement[2] << " stillFrames=" << Stats.StillFrames << " stillBoth=" << Stats.StillBothRate
         << " supportAgree=" << Stats.SupportAgreement << " runs=" << Runs.size() << " minRun=" << Stats.MinRun
         << " medianRun=" << Stats.MedianRun;
    Report = Text.str();
    return true;
}

bool BakeFootContacts(const std::string& OutputPath, std::string& Report)
{
    std::vector<uint8_t> Bytes;
    FFootContactStats Stats;
    if(!BuildFootContactsFile(Bytes, Stats, Report)) return false;
    FFootContactFile File;
    std::string Reason;
    // 写出前按加载器自校验一次，保证产物可被场景接受。
    uint64_t Hash = 0;
    if(!HashFileFnv1a64(GetSandboxResourcePath("long_motion.bin"), Hash, Reason)
        || ParseFootContacts(Bytes.data(), Bytes.size(), Stats.Frames, Hash, File, Reason) != EFootContactError::None)
    { Report = "自校验失败：" + Reason; return false; }
    if(!SaveFootContacts(OutputPath, File, Reason)) { Report = Reason; return false; }
    return true;
}

#include "SceneRender.inl"
