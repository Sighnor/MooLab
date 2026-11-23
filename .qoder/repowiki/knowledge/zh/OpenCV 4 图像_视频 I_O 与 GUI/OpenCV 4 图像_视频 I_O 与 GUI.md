---
kind: external_dependency
name: OpenCV 4（仅 core/imgproc/imgcodecs）
slug: opencv4
category: external_dependency
category_hints:
    - vendor_identity
scope:
    - '**'
source_files:
    - vcpkg.json
    - CMakeLists.txt
    - Tools/ImageDiff/ImageDiff.cpp
---

> 本文档已按 2.0（分支 `Version_2.0`）现状修正：OpenCV 不再承担窗口、图像显示或滑块。

### 身份与用途
- Intel 开源计算机视觉库。2.0 中 `find_package(OpenCV CONFIG REQUIRED COMPONENTS core imgproc imgcodecs)`，vcpkg 只启用 `jpeg`/`png`/`thread` 特性，**没有 highgui、videoio、win32ui、dnn**。
- 用途收敛为两类：
  - **资源解码**：`imread` 读取 `skybox.png`、`BRDFLut.png`、`EavgLut.png`、`resources/UI/*.png` 等，转成引擎内部 `texture` 结构（见 `Source/Legacy/texture.hpp`）。
  - **截图与对拍**：`--capture-frame` 写 PNG（`Application.cpp`），`Tools/ImageDiff/ImageDiff.cpp` 读两张 PNG 做逐像素差值统计并输出 PSNR。
- 窗口显示与参数调节改由 GLFW/ImGui 承担（见 raylib 卡片中的现状说明）。

### 集成方式
- 依赖由 vcpkg manifest 提供；`MooSandbox`、`MooLab`、`MooImageDiff` 通过 `${OpenCV_LIBS}` 链接，运行库 DLL 由 applocal 部署到 exe 同目录。
- 无 `cv::namedWindow`/`imshow`/`createTrackbar`/`waitKey` 调用；主循环不依赖任何 OpenCV GUI 事件泵。

### 稳定约束
- 只使用 `core`/`imgproc`/`imgcodecs` 三个组件，新增组件必须同步改 `CMakeLists.txt` 与 `vcpkg.json`。
- 截图与对拍的量化规则（8 位舍入、每像素最大通道差、≤3 比例 ≥99%、>32 比例 ≤0.5%）依赖 OpenCV 的解码/编码行为；升级 OpenCV 版本后需重跑 `verify_gl.py capture` 复核基准哈希（README 记录的是 4.14.0 实测值）。
