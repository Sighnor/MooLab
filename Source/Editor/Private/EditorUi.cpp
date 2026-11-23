#include "EditorUi.h"
#include <algorithm>
#include <imgui.h>
#include <imgui_internal.h>

FEditorUi::FEditorUi(bool bChineseFont) : bChinese(bChineseFont) {}
const char* FEditorUi::Label(const char* Chinese, const char* English) const { return bChinese ? Chinese : English; }
bool FEditorUi::EditMaterial(FSurfaceMaterial& Material) const
{
    // 范围与场景校验一致；越界输入由滑条钳制，场景仍会整体拒绝非法值
    bool bChanged = ImGui::SliderFloat(Label("金属度", "Metallic"), &Material.Metallic, 0.f, 1.f, "%.3f", ImGuiSliderFlags_AlwaysClamp);
    bChanged |= ImGui::SliderFloat(Label("粗糙度", "Roughness"), &Material.Roughness, 0.05f, 1.f, "%.3f", ImGuiSliderFlags_AlwaysClamp);
    return bChanged;
}
void FEditorUi::RestoreLayout(unsigned int DockId)
{
    ImGui::DockBuilderRemoveNode(DockId);
    ImGui::DockBuilderAddNode(DockId, ImGuiDockNodeFlags_DockSpace);
    ImGui::DockBuilderSetNodeSize(DockId, ImGui::GetMainViewport()->WorkSize);
    ImGuiID Center = DockId, Left = 0, Right = 0, Bottom = 0;
    ImGui::DockBuilderSplitNode(Center, ImGuiDir_Left, 0.16f, &Left, &Center);
    ImGui::DockBuilderSplitNode(Center, ImGuiDir_Right, 0.26f, &Right, &Center);
    ImGui::DockBuilderSplitNode(Center, ImGuiDir_Down, 0.25f, &Bottom, &Center);
    ImGui::DockBuilderDockWindow(Label("场景大纲", "Outliner"), Left);
    ImGui::DockBuilderDockWindow("Viewport", Center);
    ImGui::DockBuilderDockWindow(Label("属性", "Properties"), Right);
    ImGui::DockBuilderDockWindow(Label("播放 / 仿真", "Playback / Physics"), Bottom);
    ImGui::DockBuilderDockWindow("Motion Matching", Bottom);
    ImGui::DockBuilderDockWindow(Label("统计", "Statistics"), Bottom);
    ImGui::DockBuilderFinish(DockId);
    bViewportOpen = true;
    bRestoreLayout = false;
    FocusStage = 2;
}
FEditorFrame FEditorUi::Build(const FSceneSettings& Settings, const FSceneSnapshot& Snapshot,
                             const FFrameStats& Stats, const FViewportTexture& Texture)
{
    FEditorFrame Result;
    Result.Settings = Settings;
    auto& Value = Result.Settings;
    auto& Command = Result.Command;
    const bool bMatching = Snapshot.Mode == ESceneMode::MotionMatching;
    bool bChanged = false;
    if(ImGui::BeginMainMenuBar())
    {
        if(ImGui::Button(Label("物理演示", "Physics"))) { Command.bSetMode = true; Command.Mode = ESceneMode::Physics; }
        ImGui::SameLine();
        ImGui::BeginDisabled(!Snapshot.bMatchingAvailable);
        if(ImGui::Button(Label("角色控制 / Motion Matching", "Character / Motion Matching")))
        { Command.bSetMode = true; Command.Mode = ESceneMode::MotionMatching; }
        ImGui::EndDisabled();
        ImGui::SameLine();
        if(ImGui::Button(Label("恢复布局", "Restore Layout"))) bRestoreLayout = true;
        ImGui::SameLine();
        ImGui::TextUnformatted(Label("WASD行走 · 方向键相机 · 空格暂停 · .单步（先聚焦视口）",
                                     "WASD move / arrows camera / space pause / . step (focus Viewport)"));
        ImGui::EndMainMenuBar();
    }
    const ImGuiID DockId = ImGui::GetID("MooLabDock");
    if(bRestoreLayout) RestoreLayout(DockId);
    ImGui::DockSpaceOverViewport(DockId, ImGui::GetMainViewport());
    if(ImGui::Begin(Label("场景大纲", "Outliner")))
    {
        // 固定项 Camera/PointLight/Ground/Character 之后依次是全部胶囊。
        auto Item = [&](const char* Text, EEditorObject Type, int CapsuleIndex)
        {
            const bool bSelected = Selection.Type == Type && (Type != EEditorObject::Capsule || Selection.CapsuleIndex == CapsuleIndex);
            if(ImGui::Selectable(Text, bSelected)) Selection = {Type, CapsuleIndex};
        };
        Item("Camera", EEditorObject::Camera, -1);
        Item("PointLight", EEditorObject::PointLight, -1);
        Item("Ground", EEditorObject::Ground, -1);
        Item("Character", EEditorObject::Character, -1);
        for(int Index = 0; Index < static_cast<int>(Snapshot.Bones.size()); ++Index)
        {
            ImGui::PushID(Index);
            char Text[64];
            ImFormatString(Text, sizeof(Text), "Capsule %02d", Index);
            Item(Text, EEditorObject::Capsule, Index);
            ImGui::PopID();
        }
    }
    ImGui::End();
    if(bViewportOpen)
    {
        if(ImGui::Begin("Viewport", &bViewportOpen))
        {
            const ImVec2 Available = ImGui::GetContentRegionAvail();
            const bool bValidTexture = Texture.Handle != 0 && Texture.Width > 0 && Texture.Height > 0;
            const float Scale = bValidTexture
                ? std::max(0.f, std::min(Available.x / Texture.Width, Available.y / Texture.Height)) : 0.f;
            const ImVec2 Size(Texture.Width * Scale, Texture.Height * Scale);
            const ImVec2 Cursor = ImGui::GetCursorPos();
            ImGui::SetCursorPos(ImVec2(Cursor.x + (Available.x - Size.x) * 0.5f, Cursor.y + (Available.y - Size.y) * 0.5f));
            if(Scale > 0.f) ImGui::Image(static_cast<ImTextureID>(Texture.Handle), Size,
                                      ImVec2(Texture.U0, Texture.V0), ImVec2(Texture.U1, Texture.V1));
            Result.bAllowSceneInput = ImGui::IsWindowFocused(ImGuiFocusedFlags_RootAndChildWindows);
        }
        ImGui::End();
    }
    if(ImGui::Begin(Label("属性", "Properties")))
    {
        const int CapsuleCount = static_cast<int>(std::min(Snapshot.Bones.size(), Value.CapsuleMaterials.size()));
        if(Selection.Type == EEditorObject::Capsule && (Selection.CapsuleIndex < 0 || Selection.CapsuleIndex >= CapsuleCount))
            Selection = {};
        if(Selection.Type == EEditorObject::Camera)
        {
            ImGui::Text("Camera: %.3f %.3f %.3f", Snapshot.CameraPosition[0], Snapshot.CameraPosition[1], Snapshot.CameraPosition[2]);
            bChanged |= ImGui::SliderFloat(Label("跟随距离", "Distance"), &Value.CameraDistance, 0.2f, 20.f);
            bChanged |= ImGui::SliderFloat(Label("跟随高度", "Height"), &Value.CameraHeight, -5.f, 5.f);
            bChanged |= ImGui::SliderFloat(Label("转向灵敏度", "Sensitivity"), &Value.CameraSensitivity, 0.f, 0.2f);
        }
        else if(Selection.Type == EEditorObject::PointLight)
        {
            bChanged |= ImGui::Checkbox(Label("跟随角色", "Follow character"), &Value.bFollowLight);
            ImGui::BeginDisabled(Value.bFollowLight);
            bChanged |= ImGui::SliderFloat3(Label("手动位置", "Manual position"), Value.LightPosition.data(), -100.f, 100.f);
            ImGui::EndDisabled();
            bChanged |= ImGui::SliderFloat(Label("辐射基值", "Radiance base"), &Value.LightBase, 0.f, 100.f);
            bChanged |= ImGui::SliderFloat(Label("辐射振幅", "Radiance amplitude"), &Value.LightAmplitude, 0.f, Value.LightBase);
            ImGui::Text("Position: %.2f %.2f %.2f", Snapshot.LightPosition[0], Snapshot.LightPosition[1], Snapshot.LightPosition[2]);
            ImGui::Text("Radiance: %.3f", Snapshot.LightRadiance[0]);
        }
        else if(Selection.Type == EEditorObject::Ground)
        {
            bChanged |= ImGui::Checkbox(Label("显示地面", "Show ground"), &Value.Ground.bVisible);
            bChanged |= ImGui::SliderFloat(Label("格长 m", "Cell size m"), &Value.Ground.CellSize, 0.05f, 10.f, "%.3f",
                                           ImGuiSliderFlags_AlwaysClamp | ImGuiSliderFlags_Logarithmic);
            bChanged |= EditMaterial(Value.Ground.Material);
            const auto& Plane = Snapshot.GroundPlane;
            if(bMatching)
            {
                bChanged |= ImGui::SliderFloat(Label("地面高度 m", "Ground height m"), &Value.Ground.MatchingHeight, -2.f, 2.f, "%.3f",
                                               ImGuiSliderFlags_AlwaysClamp);
            }
            else
            {
                ImGui::Text(Label("实际平面 A/B/C/D：%.4f %.4f %.4f %.4f", "Active plane A/B/C/D: %.4f %.4f %.4f %.4f"), Plane[0], Plane[1], Plane[2], Plane[3]);
                ImGui::TextWrapped("%s", Label("Physics下地面贴合接触平面，请在“播放 / 仿真”面板修改。",
                                               "In Physics the ground follows the contact plane; edit it in Playback / Physics."));
            }
            ImGui::TextWrapped("%s", Label("地面格子固定在世界中，用于观察滑步；不参与仿真或匹配。",
                                           "World-fixed grid for observing foot sliding; it does not affect simulation or matching."));
        }
        else
        {
            const bool bCapsule = Selection.Type == EEditorObject::Capsule;
            bChanged |= EditMaterial(bCapsule ? Value.CapsuleMaterials[Selection.CapsuleIndex] : Value.CharacterMaterial);
            bChanged |= ImGui::Checkbox(Label("显示胶囊", "Show capsules"), &Value.bShowCapsules);
            bChanged |= ImGui::Checkbox(Label("显示接触点", "Show contacts"), &Value.bShowContacts);
            bChanged |= ImGui::Checkbox(Label("显示角色皮肤", "Show character skin"), &Value.bShowCharacter);
            if(bCapsule)
            {
                const auto& Bone = Snapshot.Bones[Selection.CapsuleIndex];
                ImGui::Text("Position: %.3f %.3f %.3f", Bone.Position[0], Bone.Position[1], Bone.Position[2]);
                ImGui::Text("Rotation: %.3f %.3f %.3f %.3f", Bone.Rotation[0], Bone.Rotation[1], Bone.Rotation[2], Bone.Rotation[3]);
                ImGui::Text("Contact: %s", Bone.bContact ? "true" : "false");
                ImGui::TextUnformatted(Label("位姿只读（由场景驱动）", "Pose is read-only (scene driven)"));
            }
            else
            {
                ImGui::Text("Root: %.3f %.3f %.3f", Snapshot.RootPosition[0], Snapshot.RootPosition[1], Snapshot.RootPosition[2]);
                ImGui::Text("Yaw: %.3f / Speed: %.3f m/s", Snapshot.RootYaw, Snapshot.RootSpeed);
            }
        }
        if(!Snapshot.Message.empty()) ImGui::TextWrapped("%s", bChinese ? Snapshot.Message.c_str() : "Command rejected or Motion Matching unavailable; see scene diagnostics.");
    }
    ImGui::End();
    // 布局恢复后先选中底部“播放”页签，下一帧把焦点交还 Viewport。
    if(FocusStage == 2) ImGui::SetNextWindowFocus();
    else if(FocusStage == 1) ImGui::SetWindowFocus("Viewport");
    if(FocusStage > 0) --FocusStage;
    if(ImGui::Begin(Label("播放 / 仿真", "Playback / Physics")))
    {
        if(ImGui::Button(Snapshot.bPaused ? Label("继续", "Resume") : Label("暂停", "Pause"))) Command.bTogglePause = true;
        ImGui::SameLine(); if(ImGui::Button(Label("单步", "Step"))) Command.bSingleStep = true;
        ImGui::SameLine(); if(ImGui::Button(Label("重置", "Reset"))) Command.bReset = true;
        ImGui::SameLine(); if(ImGui::Button(Label("恢复默认", "Defaults"))) Command.bRestoreDefaults = true;
        ImGui::BeginDisabled(bMatching);
        bChanged |= ImGui::InputDouble(Label("重力 0..30", "Gravity 0..30"), &Value.Gravity, 0.1, 1., "%.4f");
        bChanged |= ImGui::InputDouble(Label("阻尼 0..5", "Damping 0..5"), &Value.Damping, 0.01, 0.1, "%.4f");
        bChanged |= ImGui::SliderInt(Label("子步", "Substeps"), &Value.Substeps, 1, 8);
        int Control = static_cast<int>(Value.Control);
        if(ImGui::Combo(Label("力矩模式", "Control"), &Control, "Global\0Local\0Off\0"))
        { Value.Control = static_cast<EControlMode>(Control); bChanged = true; }
        bChanged |= ImGui::InputInt(Label("控制截止步", "Control end frame"), &Value.ControlEndFrame);
        bChanged |= ImGui::InputScalarN(Label("接触平面 A/B/C/D", "Plane A/B/C/D"), ImGuiDataType_Double, Value.ContactPlane.data(), 4, nullptr, nullptr, "%.4f");
        bChanged |= ImGui::InputDouble(Label("接触容差", "Contact epsilon"), &Value.ContactEpsilon, 0.001, 0.01, "%.4f");
        ImGui::EndDisabled();
        if(!bMatching) ImGui::TextUnformatted(Label("Physics下WASD不驱动求解；请启用角色控制。", "WASD does not drive Physics; enable Character control."));
    }
    ImGui::End();
    if(ImGui::Begin("Motion Matching"))
    {
        ImGui::Text("Available: %s / Candidates: %d / Excluded: %d", Snapshot.bMatchingAvailable ? "yes" : "no",
                    Snapshot.CandidateCount, Snapshot.ExcludedCount);
        ImGui::BeginDisabled(!bMatching);
        bChanged |= ImGui::SliderFloat(Label("速度上限 m/s", "Max speed m/s"), &Value.MoveSpeed, 0.f, 10.f);
        bChanged |= ImGui::SliderFloat3(Label("速度 / 转向 / 脚 权重", "Velocity / yaw / feet weights"), Value.MatchingWeights.data(), 0.f, 10.f);
        bChanged |= ImGui::Checkbox(Label("脚锁IK", "Foot lock IK"), &Value.bFootIk);
        ImGui::EndDisabled();
        for(int Side = 0; Side < 2; ++Side)
        {
            const FFootView& Foot = Snapshot.Feet[Side];
            ImGui::Text("%s: %s / %s / %s %.3f m%s", Side == 0 ? Label("左脚", "Left foot") : Label("右脚", "Right foot"),
                        Foot.bContactLabel ? Label("触地", "contact") : Label("悬空", "air"),
                        Foot.bLocked ? Label("锁定", "locked") : (Foot.bFading ? Label("淡出", "fading") : Label("自由", "free")),
                        Label("残差", "residual"), Foot.ClampError,
                        Foot.bConflict ? Label(" / 膝角优先", " / knee priority") : "");
        }
        ImGui::Text("Source: %d / Clip: %d%s / Last search step: %d / Searches: %d", Snapshot.SourceFrame, Snapshot.ClipStart,
                    Snapshot.bBlending ? " (blend)" : "", Snapshot.SearchFrame, Snapshot.SearchCount);
        ImGui::Text("Cost: %.6f / Search: %.3f ms / Foot IK: %.3f ms", Snapshot.SearchCost, Snapshot.SearchMs, Snapshot.FootIkMs);
        if(!bMatching) ImGui::TextUnformatted(Label("当前Physics模式：匹配控制未启用。", "Physics mode: matching controls inactive."));
    }
    ImGui::End();
    if(ImGui::Begin(Label("统计", "Statistics")))
    {
        ImGui::Text("Cost: %.3f ms / GPU: %.3f ms / Interval: %.3f ms", Stats.CostMs, Stats.GpuMs, Stats.IntervalMs);
        ImGui::Text("Mode: %s / App: %llu / Simulation: %d", bMatching ? "Motion Matching" : "Physics",
                    static_cast<unsigned long long>(Stats.AppFrame), Snapshot.NextFrame);
        if(!bChinese) ImGui::TextUnformatted("Chinese system font unavailable: using English labels.");
    }
    ImGui::End();
    Result.bSettingsChanged = bChanged;
    Result.bAllowSceneInput &= !ImGui::GetIO().WantTextInput && !ImGui::IsAnyItemActive()
        && !ImGui::IsPopupOpen(nullptr, ImGuiPopupFlags_AnyPopupId | ImGuiPopupFlags_AnyPopupLevel);
    return Result;
}
