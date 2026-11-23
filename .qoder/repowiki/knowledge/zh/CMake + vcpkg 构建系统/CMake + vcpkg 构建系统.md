---
kind: build_system
name: CMake + vcpkg 构建系统
category: build_system
scope:
    - '**'
source_files:
    - CMakeLists.txt
    - CMakePresets.json
    - vcpkg.json
---

> 本文档已按 2.0（分支 `Version_2.0`）现状修正：单一 `main.cpp` + 头文件库的 1.x 结构已不存在，仓库无 Makefile、Dockerfile 或 CI 脚本。

## 1. 构建系统
`CMake 3.21+`（`CMakeLists.txt` 声明的最低版本）+ `CMakePresets.json`（schema v6）+ vcpkg manifest。Windows-only，preset 条件为 `hostSystemName == "Windows"`。

## 2. 关键文件
- `CMakeLists.txt`：全局编译选项、`find_package`、模块目标、着色器与工具的注册、公共头文件名检查。
- `CMakePresets.json`：`windows-base` 提供生成器/架构/工具链，`vs2026`（`Visual Studio 18 2026`）与 `vs2022`（`Visual Studio 17 2022`）为配置预设，`*-debug`/`*-release` 为构建预设；`binaryDir` 为 `build/${presetName}`。
- `vcpkg.json`：`version-string` 2.0 + `builtin-baseline` + 4 个 port（见 vcpkg 依赖管理卡片）。
- `build/`：生成目录（整体被 `.gitignore` 排除）。

## 3. 依赖解析（第 8–12 行）
```
find_package(OpenCV CONFIG REQUIRED COMPONENTS core imgproc imgcodecs)
find_package(OpenMP REQUIRED COMPONENTS CXX)
find_package(glfw3 CONFIG REQUIRED)
find_package(glad CONFIG REQUIRED)
find_package(imgui CONFIG REQUIRED)
```
全部走 CONFIG 模式，不手写 include/lib 路径；OpenMP 用 MSVC 内置实现，不额外引入 vcpkg 包。

## 4. 目标与模块
- `MooBuildOptions`（INTERFACE）：MSVC 下定义 `NOMINMAX`、`WIN32_LEAN_AND_MEAN`、`_CRT_SECURE_NO_WARNINGS`，编译选项 `/utf-8 /permissive- /W3 /MP`。
- `MooCore`（INTERFACE）：只导出 `Source/Runtime/Core/Public`，不含 GL/窗口/旧数学头。
- `moo_module(Target Directory ...)`：`add_library(... STATIC)` + 导出 `<Directory>/Public` + 链接 `MooCore` + `source_group`。
- 静态库：`MooRHI`（链接 `glad::glad`，定义 `MOOLAB_SHADER_DIR`）、`MooRenderer`、`MooApplicationCore`（`glfw`/`glad::glad`/`imgui::imgui`，Win32 下再链 `Xinput winmm`）、`MooAnimation`（纯数据，不依赖 Legacy 与 GL）、`MooSandbox`（`PRIVATE Source/Legacy`、`${OpenCV_LIBS}`、`OpenMP::OpenMP_CXX`、`MOOLAB_RESOURCE_DIR`）、`MooEditor`（只依赖 Sandbox 的 Public 与 ImGui）。
- 可执行：`MooLab`（`Source/Launch/Private/Main.cpp` + `Application.*`），链接 `MooApplicationCore MooSandbox MooRenderer MooEditor ${OpenCV_LIBS}`，并把 `Shaders/*` 注册为源文件（`source_group("Shaders")`），设置 `VS_STARTUP_PROJECT` 与 `VS_DEBUGGER_WORKING_DIRECTORY=resources`。
- 工具：`MooImageDiff` 参与默认构建；`SceneProbe`、`FootContactBake`、`gl_probe`、`MaterialProbe` 为 `EXCLUDE_FROM_ALL`，需显式 `--target` 构建。

## 5. 约定与约束
- **公共头文件名全仓唯一**：`MOO_PUBLIC_HEADERS` 列表在配置期逐个取 basename，重复即 `message(FATAL_ERROR "Duplicate public header basename")`；新增公共头必须同时登记到该列表。
- **新增文件/目录需显式注册**：`moo_module` 不做 glob；新增模块要写 `moo_module(...)` 并补导出路径。
- **Legacy 隔离**：`Source/Legacy/*.hpp` 含大量非 inline 定义，只允许 `Source/Sandbox/Private/SandboxScene.cpp` 包含；Sandbox 下的 `.inl` 标为 `HEADER_FILE_ONLY`，仅由该 .cpp 包含，避免 ODR 重复定义。
- **编译期路径宏**：`MOOLAB_RESOURCE_DIR`（MooSandbox）与 `MOOLAB_SHADER_DIR`（MooRHI）都是源码目录的绝对路径，与启动工作目录无关。
- **字符编码**：源码为无 BOM UTF-8（含中文注释），MSVC 必须保留 `/utf-8`。
- **标准**：`CMAKE_CXX_STANDARD 17`、`STANDARD_REQUIRED ON`、`EXTENSIONS OFF`，并开 `CMAKE_EXPORT_COMPILE_COMMANDS`。
- **生成器**：VS 2026 需用 VS 自带的 CMake（4.3+），系统盘上的 CMake 3.27 不认识 `Visual Studio 18 2026`。
- **版本号**：仅 `vcpkg.json` 的 `version-string` 声明，CMake 中不重复。
