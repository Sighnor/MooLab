# MooLab 2.0

角色运动渲染与仿真器  
依赖 GLFW3、glad、ImGui、OpenGL 4.6、OpenCV（仅 core/imgproc/imgcodecs）与 OpenMP；系统库 XInput 负责手柄输入，winmm 仅用于计时器回退。raylib 与 OpenCV highgui 不再用于运行。
主要参考了GAMES101, GAMES202, GAMES105, https://github.com/orangeduck/Motion-Matching, https://github.com/raysan5/raylib  
用的素材也来自上述项目，感谢各位大佬  

> **版本**：2.0（分支 `Version_2.0`）。相对 1.x（提交 `fadced7`，v1.09）的主要变化是新增默认 GPU 渲染路径、停靠式编辑器与 `Source/` 模块化目录重组；1.x 不再维护，历史代码可在仓库旧提交中查看。逐条改动见 [CHANGELOG.md](CHANGELOG.md)。

## 构建（Windows / Visual Studio）
依赖通过 vcpkg manifest（`vcpkg.json`）自动安装，无需手动配置 OpenCV/GLFW/glad/ImGui；运行需要支持 OpenGL 4.6 的驱动（CPU 参考模式也用 GL 窗口呈现）。
1. 安装 Visual Studio 2026 或 2022，勾选“使用 C++ 的桌面开发”（含 CMake 与 vcpkg 组件）
2. 设置环境变量 `VCPKG_ROOT`，可直接指向 VS 自带的 vcpkg，例如 `C:\Program Files\Microsoft Visual Studio\18\Community\VC\vcpkg`
3. 任选一种方式：
   - VS 中“打开本地文件夹”选择仓库根目录，配置预设选 `vs2026`（或 `vs2022`），启动项选 `MooLab.exe`
   - 命令行：`cmake --preset vs2026` 后 `cmake --build --preset vs2026-release`，或打开 `build/vs2026/MooLab.slnx`（vs2022 预设为 `build/vs2022/MooLab.sln`）
4. 首次配置会编译依赖，耗时数分钟；产物位于 `build/<preset>/<Config>/`，运行库 DLL 自动拷贝到同目录
5. 资源和 `Shaders/` 路径由编译期宏固定到源码目录；资源目录不存在时警告并退回工作目录相对的 `../resources/`（进程工作目录不被修改），着色器打开/编译失败时明确退出（移动工程后请重新配置、编译）。性能相关测试请使用 Release；Debug 用 `cmake --build --preset vs2026-debug`。

## 目录
模块内划分 Public/Private，各模块只导出 Public 包含路径，公共头文件名全仓唯一（CMake 配置期检查）：
```text
Source/
  Runtime/Core/Public/            Core.h、GlContext.h、ViewportTexture.h（INTERFACE，无GL/窗口/旧数学头）
  Runtime/RHI/                    GlRhi（GL对象RAII、时间戳、着色器读取）
  Runtime/Renderer/               RenderTypes、GpuRenderer、CpuImagePresenter
  Runtime/ApplicationCore/        Platform（XInput/帧节流）、AppWindow（GLFW/ImGui）、Input（键鼠采集与路由）
  Sandbox/                        SceneTypes、SandboxScene（PImpl门面；Private 下 SceneRender.inl / MotionMatching.inl）
  Editor/                         EditorUi（停靠式编辑器）
  Launch/Private/                 Main.cpp、Application（窗口/帧节拍/输入路由/截图组装）
  Legacy/                         原数学/动画/仿真/CPU渲染 .hpp
Shaders/                          GLSL
Tools/ImageDiff/                  MooImageDiff
Tools/Verify/                     verify_gl.py、run_ui_case.py、gl_probe.cpp、SceneProbe.cpp、compile_gl_probe.cmd
resources/                        运行数据（character.bin、database.bin、long_motion.bin、skybox.png、foot_contacts.bin 等必需文件）
.qoder/                           Qoder 智能体定义（agents）与 repowiki 知识库
```
- Legacy 头含大量非 inline 定义，只允许被 `Source/Sandbox/Private/SandboxScene.cpp` 包含（`.inl` 仅由它包含，CMake 标记 `HEADER_FILE_ONLY`）；其他编译单元包含 Legacy 会 ODR 多重定义。
- 依赖方向：Renderer → RHI/Core；ApplicationCore → Core；Sandbox → Renderer/Core；Editor 只用 Sandbox 公共纯数据与 ImGui；Launch 组装生成 `MooLab.exe`。

## 运行
- 窗口为 1600×900 停靠式编辑器：左侧场景大纲（Camera、PointLight、Character、22 个胶囊），中间 Viewport（离屏 720×720 等比留边显示），右侧属性，底部播放/仿真、Motion Matching、统计页签；可“恢复布局”，不写 imgui.ini。中文标签使用系统已装字体，加载失败退回英文。
- 场景模式：默认 `Physics`（原刚体+力矩控制演示，图像与历史一致）；工具栏可切到“角色控制 / Motion Matching”。MM 只读 `resources/long_motion.bin` 与 `database.bin` 在内存中搜索与混合，不改写资源；数据校验失败时禁止进入 MM，Physics 仍可用。
- 键盘：WASD 按相机水平朝向控制角色（S 为转身向后走，不是倒走），方向键环绕相机；空格暂停/继续，`.` 单步。仅在窗口与 Viewport 聚焦、且没有输入框占用时生效；鼠标只操作 UI。键盘与 XInput 手柄（左摇杆移动、右摇杆相机）贡献先相加再做径向死区。Physics 模式下左摇杆不改变物理求解。
- 可调参数：重力、阻尼、子步、接触平面与 eps、力矩模式与控制截止步、相机跟随距离/高度/灵敏度、灯光跟随/位置/辐射、胶囊 roughness、显示开关、MM 速度上限与匹配权重；非法值拒绝并显示原因。“恢复默认”再“重置”可还原历史初态。
- ESC 或关闭窗口退出；目标帧率 60，优先用高精度可等待计时器节流，系统不支持时回退到 1ms 计时器精度的 sleep
- 默认 GPU 路径：光源立方体、22 个 PBR 胶囊与绿色接触点，顶点矩阵变换/光栅/片元着色在 GPU 执行；FK/蒙皮/仿真仍在 CPU，角色皮肤默认不绘制（属性中可开启）。
- 控制台输出 `Cost`（窗口消息、编辑器 UI、仿真、渲染、显示提交及命中截图写盘，不含节流等待）、`Interval`（上一帧完整间隔）与 `Gpu`（仅离屏渲染的 GPU 时间戳，不含呈现/UI/读回，滞后若干帧，未就绪为 0；CPU 模式恒 0）。GPU 提交异步，不能用 Cost 当作 GPU 执行时间。退出时输出平均/最大帧间隔；日志整行写出，不逐帧强制刷新。
- 命令行参数（用于自动验证）：
  - `--frames N`：运行 N 帧后退出；缺少值时非零退出
  - `--capture-frame N PATH`：第 N 帧保存为图片（从 0 开始），可重复，不含 ImGui；参数不完整、写入失败/异常或未达到请求帧均非零退出
  - `--renderer gpu|cpu`：默认 gpu，cpu 使用原软光栅与截图量化规则；非法模式非零退出
  - `--scene physics|motion-matching`：初始场景模式，默认 physics；非法值或缺值非零退出
  - `--hidden`：创建隐藏 GLFW 窗口用于无人值守对拍，仍需要 GL 上下文，不是无显示服务的 headless
  - `--unlimited-fps`：取消 60 帧限制
- 对拍例（仓库根目录 PowerShell）：
  ```powershell
  & .\build\vs2026\Release\MooLab.exe --renderer cpu --hidden --frames 301 --unlimited-fps --capture-frame 60 cpu_060.png --capture-frame 300 cpu_300.png
  & .\build\vs2026\Release\MooLab.exe --renderer gpu --hidden --frames 301 --unlimited-fps --capture-frame 60 gpu_060.png --capture-frame 300 gpu_300.png
  & .\build\vs2026\Release\MooImageDiff.exe gpu_060.png cpu_060.png
  & .\build\vs2026\Release\MooImageDiff.exe gpu_300.png cpu_300.png
  ```
  MooImageDiff 默认以每像素最大通道差统计：差 ≤3 的比例 ≥99%，差 >32 的比例 ≤0.5%；退出码 0=通过、1=阈值失败、2=输入错误、3=尺寸不一致。可用 `--diff-out diff.png` 输出差异图。负对照 `gpu_060.png` 对 `cpu_300.png` 必须失败。未接手柄时 Release 的 CPU 60/300 帧 SHA-256 前16位为 `F32E92295DFD7FFF` / `6CC141F39195243E`（本机 OpenCV 4.14.0 实测；`foot_contacts.bin` 重新烘焙后姿态已变，与 1.x 的 `004F883938E6AD0A` / `EB3E20DCBC32FF56` 不同）；CPU 单/多线程、GPU 重跑与 hidden/可见同帧截图应逐字节一致。
- 自动验证（需 Python 与 Pillow，输出在 `build/vs2026`，可加 `--build-dir=...`）：`python Tools\Verify\verify_gl.py capture`（7 组截图、历史哈希、对拍与 CLI 退出码）、`probe`（GL 探针与 4 个故障用例，先运行 `Tools\Verify\compile_gl_probe.cmd`）、`perf`（3×600 帧性能）。场景探针 `cmake --build --preset vs2026-release --target SceneProbe` 后运行 `build\vs2026\Release\SceneProbe.exe <输出前缀>`，检查 Physics 重置/模式往返的历史帧、暂停/单步、MM 边界与 30000 步长测及运动预算（0=通过、1=硬失败、3=预算失败）。
- GPU 保留 CPU 的非透视插值、纹理 UV 互换/最近邻、F0=albedo 与观察空间深度等非标准行为；对拍通过证明迁移等价，不证明 PBR 物理正确。
- 当前场景的低CPU启动方式（PowerShell，策略为可选，程序未写死默认）：
  ```powershell
  $env:OMP_WAIT_POLICY='PASSIVE'
  & .\build\vs2026\Release\MooLab.exe
  ```
  本机 i9-12900KF / RTX 5080（NVIDIA 610.88）/ Release / 720×720 / `--hidden` / 60 帧节流 / `OMP_WAIT_POLICY=PASSIVE` / 默认线程数，每种模式 3×600 帧（`Cost`、`Gpu` 排除前 10 帧）：
  - 三次均值（编辑器版本，Physics 场景，`Cost` 含编辑器 UI 构建）：GPU 路径 `Cost` P50/P99 约 2.91/4.08ms，离屏 `Gpu` P50 约 0.014ms；CPU 路径 `Cost` P50/P99 约 8.40/10.15ms（CPU 模式仅在截图帧做 8 位转换）。
  - 进程 CPU 秒/墙钟秒（含加载启动）GPU 约 0.25、CPU 约 1.01 个等效核；两种模式完整帧间隔都约 17ms，受 60 帧节流限制，不能推出 FPS 提升倍数。
  - Motion Matching 由 `SceneProbe` 单线程测得平均每次全量搜索约 0.17ms、平均每步约 0.18ms（Release，30000 步）；未纳入上述窗口性能口径。
  - 场景仅为光源立方体、22 个胶囊与接触点；未测不限帧率吞吐、可见窗口性能及分段热点。等效核不是物理核数，不能外推为其他场景的保证。

  旧 OpenCV highgui 窗口时期同场景 CPU 路径曾测得约 1.15~1.21 个等效核；未设置等待策略时一次实测约 18.8 个等效核。该历史数据窗口层不同，不能与上面的 GL 数据直接横比。可用 `Remove-Item Env:\OMP_WAIT_POLICY` 恢复运行库默认策略。

效果展示：
使用全局力矩控制，一段时间后撤去力矩(gif压制帧率有点低)
![1](https://github.com/Sighnor/MooLab/assets/109077418/7d2fed0c-a555-40aa-8ca1-fa1b95a19b28)


