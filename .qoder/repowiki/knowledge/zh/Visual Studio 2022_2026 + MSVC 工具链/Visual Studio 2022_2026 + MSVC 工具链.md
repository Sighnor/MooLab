---
kind: external_dependency
name: Visual Studio 2022/2026 + MSVC 工具链
slug: visual-studio-msvc
category: external_dependency
category_hints:
    - client_constraint
scope:
    - '**'
source_files:
    - CMakePresets.json
    - CMakeLists.txt
---

### 身份与用途
- 官方支持的唯一开发环境：VS 2022（`Visual Studio 17 2022`）或 VS 2026（`Visual Studio 18 2026`），对应 preset `vs2022` / `vs2026`。
- 项目历史上从 MinGW/VSCode 迁到 MSVC，不再支持 GCC/Clang 工具链，也没有 Linux/macOS 分支。

### 工具链注意事项
- **CMake 版本要与生成器匹配**：`Visual Studio 18 2026` 需要 VS 自带的 CMake（4.3.x，如 `<VS>\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe`）；系统单独安装的 CMake 3.27 不认识该生成器，会报生成器无效。
- **vcpkg 用 VS 自带的即可**：`<VS安装目录>\VC\vcpkg`，把 `VCPKG_ROOT` 指过去（可写 User 级环境变量）。
- **运行时需要 OpenGL 4.6 驱动**：默认 GPU 路径与 CPU 参考路径都通过 GL 窗口呈现，`--hidden` 只是隐藏窗口，不是无显示服务。

### 稳定约束
- 仅 Windows + x64；triplet 固定 `x64-windows`，architecture 固定 `x64`。
- MSVC 编译选项 `/utf-8 /permissive- /W3 /MP`，宏 `NOMINMAX`、`WIN32_LEAN_AND_MEAN`、`_CRT_SECURE_NO_WARNINGS`；源码保持无 BOM UTF-8（含中文注释与中文 UI 文本）。
- OpenMP 用 MSVC 内置实现（`find_package(OpenMP CONFIG … components CXX)`），不引入独立 vcpkg 包。
- 源码对 MSVC 做过最小兼容改动：`__restrict__` → `MOO_RESTRICT` 宏，变长数组 → `std::vector`。
- 系统库依赖：`Xinput`（手柄）、`winmm`（计时器精度回退），随 `MooApplicationCore` 在 Win32 下链接。
