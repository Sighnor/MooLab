---
kind: external_dependency
name: vcpkg 依赖管理
slug: vcpkg
category: external_dependency
category_hints:
    - vendor_identity
scope:
    - '**'
source_files:
    - vcpkg.json
    - CMakePresets.json
---

> 本文档已按 2.0（分支 `Version_2.0`）现状修正，1.x 的 raylib / OpenCV highgui 描述作废。

### 身份与用途
- Microsoft 的 C/C++ 包管理器，是 2.0 唯一的第三方依赖来源（无 git submodule、无手工放入的第三方库）。
- `vcpkg.json` 声明 4 个 port：`opencv4`（仅 `core`/`imgproc`/`imgcodecs`，特性 `jpeg`/`png`/`thread`）、`glfw3`、`glad`（`gl-api-46`）、`imgui`（`glfw-binding`/`opengl3-binding`/`docking-experimental`）。
- raylib 与 OpenCV highgui 已从依赖中移除：窗口与上下文改用 GLFW + glad，UI 改用 ImGui，手柄改用系统 XInput，图像呈现改用 GL 帧缓冲 + ImGui Viewport。

### 集成方式
- `CMakePresets.json` 的 `windows-base` 把 `toolchainFile` 指向 `$env{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake`，并固定 `VCPKG_TARGET_TRIPLET=x64-windows`。
- 环境变量 `VCPKG_ROOT` 必须指向 vcpkg 根目录，可直接用 VS 自带的 `<VS安装目录>\VC\vcpkg`；当前开发机已写入 User 级环境变量（新开 shell 生效）。
- 依赖安装到仓库下的 `vcpkg_installed/`（已被 `.gitignore` 排除），运行库 DLL 由 vcpkg 的 applocal 部署机制自动拷贝到 exe 同目录。

### 稳定约束
- `builtin-baseline` 固定为 `b3df23b9d6fd69a3ce8fe8ccf943079bcc1499ad`，不随上游默认分支漂移。
- `version-string`（2.0）只作标识，不参与依赖 ABI 计算；实测只改版本号再配置会输出 "All requested installations completed successfully"，不会重装依赖。
- 仅 Windows + MSVC + x64 一种组合，无 Linux/macOS/ARM64 分支。
- 首次配置需联网拉取 registry 与源码；在受限网络下需给 git 配代理（`http.proxy`），大 pack 经代理传输挂起时追加 `http.version=HTTP/1.1` 与 `http.lowSpeedLimit`。
