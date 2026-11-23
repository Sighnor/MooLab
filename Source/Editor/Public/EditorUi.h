#pragma once
#include "SceneTypes.h"
#include "ViewportTexture.h"
#include <cstdint>

struct FEditorFrame
{
    FSceneSettings Settings;
    FSceneCommand Command;
    bool bSettingsChanged = false;
    bool bAllowSceneInput = false;
};
class FEditorUi
{
public:
    explicit FEditorUi(bool bChineseFont);
    FEditorFrame Build(const FSceneSettings& Settings, const FSceneSnapshot& Snapshot,
                       const FFrameStats& Stats, const FViewportTexture& Texture);
private:
    // 选择按对象类型记录，与大纲显示位置脱钩；CapsuleIndex仅对Capsule有效，访问前须检查范围。
    enum class EEditorObject : uint8_t { Camera, PointLight, Ground, Character, Capsule };
    struct FEditorSelection
    {
        EEditorObject Type = EEditorObject::Character;
        int CapsuleIndex = -1;
    };
    const char* Label(const char* Chinese, const char* English) const;
    void RestoreLayout(unsigned int DockId);
    bool EditMaterial(FSurfaceMaterial& Material) const;
    bool bChinese = false;
    bool bRestoreLayout = true;
    bool bViewportOpen = true;
    int FocusStage = 0;
    FEditorSelection Selection;
};
