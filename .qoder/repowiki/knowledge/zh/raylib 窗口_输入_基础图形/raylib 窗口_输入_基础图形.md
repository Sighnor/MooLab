---
kind: external_dependency
name: raylib（2.0 已移除）与 GLFW/ImGui 窗口输入层
slug: raylib
category: external_dependency
category_hints:
    - vendor_identity
scope:
    - '**'
source_files:
    - vcpkg.json
    - Source/Runtime/ApplicationCore/Public/AppWindow.h
    - Source/Runtime/ApplicationCore/Public/Platform.h
---

### 现状（2.0）
- **raylib 不再是依赖**：`vcpkg.json` 中没有该 port，代码里也不包含 raylib 头文件；README 保留致敬链接，仅作历史出处说明。
- raylib 原本承担的三件事已拆分为：
  - **窗口 + GL 上下文 + UI**：GLFW3 + glad（OpenGL 4.6 core）+ Dear ImGui（docking），见 `AppWindow.h/cpp`。
  - **手柄输入**：系统库 XInput，封装在 `Platform.cpp`（`PollGamepad`），保留原 raylib 的每轴 |v|≤0.1 归零与 Y 轴向下为正的约定；径向死区与响应曲线由调用方处理。
  - **帧率限制**：`FFrameLimiter`（`Platform.h`），优先高精度可等待计时器，系统不支持时退回 `timeBeginPeriod(1)`。

### 与窗口层相关的约束
- 每帧顺序固定为 `BeginFrame`（处理消息、开始 ImGui 帧）→ ImGui 控件与场景绘制 → 呈现到默认帧缓冲 → `EndFrame`（绘制 ImGui、交换缓冲）。
- 所有 GL 资源（`GpuRenderer`、`CpuImagePresenter`）必须在窗口之后构造、之前析构；构造时校验窗口上下文是否为当前线程的当前上下文。
- `AppWindow.h` 只前置声明 `GLFWwindow`，不暴露 glad/GLFW，可被非 GL 模块安全包含。
- 关窗或 ESC 退出；`FAppWindowDesc` 的 `bHidden` 供无人值守截图（仍需真实 GL 上下文，不是 headless）。

### 历史（1.x，已作废）
- 1.x 通过 vcpkg `raylib` port（关闭 default-features）提供窗口与 `GetGamepadAxisMovement`，并用 `SetTargetFPS(60)` 限帧；raylib 4.2→6.2 的升级记录仅作背景。
