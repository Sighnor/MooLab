---
kind: dependency_management
name: vcpkg + CMake 依赖管理
category: dependency_management
scope:
    - '**'
---

> 本文档已按 2.0（分支 `Version_2.0`）现状修正。1.x 的 raylib / OpenCV highgui 描述已全部作废。

## 1. 使用的系统/方法

本项目使用 **Microsoft vcpkg**（manifest 模式）作为 C/C++ 第三方库的包管理器，并通过 **CMake Presets** 把 vcpkg 工具链注入构建流程。2.0 的外部依赖为 `opencv4`（仅 core/imgproc/imgcodecs）、`glfw3`、`glad`（`gl-api-46`）与 `imgui`（`glfw-binding`、`opengl3-binding`、`docking-experimental`）；OpenMP、XInput、winmm 为系统级依赖，不走 vcpkg。

## 2. 关键文件

- `vcpkg.json` — vcpkg manifest：`version-string` 为 `2.0`，`builtin-baseline` 固定为 `b3df23b9d6fd69a3ce8fe8ccf943079bcc1499ad`。
- `CMakeLists.txt` — 主脚本，以 CONFIG 模式 `find_package` 消费 vcpkg 提供的 OpenCV / glfw3 / glad / imgui。
- `CMakePresets.json` — 注入 `$env{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake`，固定 `x64-windows` triplet。
- `.gitignore` — 排除 `build/`、`vcpkg_installed/`、`*.sln(x)`、`*.vcxproj` 等产物。
- `.gitattributes` — 标记 `*.bin`/`*.png`/`*.dll` 等为 binary，文本随平台。

## 3. 架构与约定

### 3.1 manifest（`vcpkg.json`）
- schema 指向 `microsoft/vcpkg-tool` 的 `vcpkg.schema.json`。
- `builtin-baseline` 锁定具体 commit，保证依赖版本可复现。
- `opencv4` 设置 `default-features: false`，仅启用 `jpeg`、`png`、`thread`；**不含 `highgui`/`win32ui`**。
- `glad` 启用 `gl-api-46`；`imgui` 启用 `glfw-binding`、`opengl3-binding`、`docking-experimental`。
- `version-string` 只作为项目元数据（2.0），不参与依赖 ABI 计算：改版本号后 reconfigure 输出 "All requested installations completed successfully"，不会触发重装。

### 3.2 CMake 集成（`CMakeLists.txt`）
- 第 8–12 行：`find_package(OpenCV CONFIG REQUIRED COMPONENTS core imgproc imgcodecs)`、`find_package(OpenMP REQUIRED COMPONENTS CXX)`、`find_package(glfw3 CONFIG REQUIRED)`、`find_package(glad CONFIG REQUIRED)`、`find_package(imgui CONFIG REQUIRED)`。
- `glad::glad`、`glfw`、`imgui::imgui` 由 `MooRHI` / `MooApplicationCore` 私有链接；`MooApplicationCore` 在 `WIN32` 下额外链接 `Xinput` 与 `winmm`。
- 不手写 `pkg-config` 或 include/lib 路径，依赖全部来自 vcpkg 安装树。
- 运行库 DLL（glfw3、opencv_core4/imgproc4/imgcodecs4、jpeg、libpng 等）由 vcpkg 的 applocal 部署自动拷贝到 `build/<preset>/<Config>/`，无需手工分发。

### 3.3 CMake Presets（`CMakePresets.json`）
- `windows-base`（hidden）：`binaryDir = ${sourceDir}/build/${presetName}`、`toolchainFile = $env{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake`、`VCPKG_TARGET_TRIPLET=x64-windows`、architecture `x64`（`strategy: set`），条件为 Windows。
- configure preset：`vs2026`（`Visual Studio 18 2026`）、`vs2022`（`Visual Studio 17 2022`）。
- build preset：`vs2026-release` / `vs2026-relwithdebinfo` / `vs2026-debug` / `vs2022-release` / `vs2022-debug`。
- `cmakeMinimumRequired` 为 3.25；`CMakeLists.txt` 自身要求 3.21。

## 4. 约定与约束

- **依赖来源单一**：所有 C/C++ 第三方库经 vcpkg 获取，不在仓库内 vendoring；`build/` 与 `vcpkg_installed/` 不入库。
- **版本锁定**：依赖版本由 `builtin-baseline` 固定，不随上游默认分支漂移。
- **最小化特性**：`opencv4` 关闭默认特性后只保留 jpeg/png/thread，显著缩短首次构建时间。
- **平台绑定**：preset 条件限定 `hostSystemName == "Windows"`，架构固定 x64，无 Linux/macOS/ARM64 配置。
- **VCPKG_ROOT**：必须指向含 `scripts/buildsystems/vcpkg.cmake` 的 vcpkg 根目录（可用 VS 自带 `…\18\Community\VC\vcpkg`）；未设置时预设的 `toolchainFile` 解析失败。
- **OpenMP 例外**：由系统 `find_package(OpenMP REQUIRED COMPONENTS CXX)` 提供（MSVC `/openmp`），不属 vcpkg 包。
- **无私有注册表**：仓库内无 `vcpkg-configuration.json` / overlay ports，依赖全部来自官方端口树。
