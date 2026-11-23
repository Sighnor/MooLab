# MooLab 变更记录

当前版本：2.0（分支 `Version_2.0`）。以提交 `fadced7`（v1.09）为对照基准，以下改动随 Version_2.0 的首次提交一并发布。

分两部分：
- 第一部分是 1.x 时期完成、此前一直未提交的 OpenGL 迁移改动。
- 第二部分是本轮"编辑器与运动控制"的改动。

---

## 一、OpenGL 迁移（1.x 时期完成）

### 构建与依赖
- 新增 `vcpkg.json`（manifest 模式）和 `CMakePresets.json`（预设 `vs2026` / `vs2022`，带 Release/Debug 构建预设）。
- 依赖换成 GLFW3、glad、ImGui、OpenGL 4.6、OpenCV（只用 core/imgproc/imgcodecs）和 OpenMP；手柄用 XInput。
- raylib 和 OpenCV highgui 不再参与运行。
- 新增 `.gitignore`，忽略 `build/` 等产物。

### 渲染
- 新增 GPU 路径并设为默认：
  - 离屏 720×720 渲染光源立方体、22 个 PBR 胶囊和接触点。
  - 顶点变换、光栅化、片元着色都在 GPU 上做；FK、蒙皮和仿真仍在 CPU。
- 保留原 CPU 软光栅作为参考模式（`--renderer cpu`），也通过 GL 窗口显示。
- GPU 路径有意复刻了 CPU 的几处非标准行为，以便逐像素对拍：
  - 非透视插值；
  - 纹理 UV 互换、最近邻采样；
  - F0 = albedo；
  - 使用观察空间深度。
- CPU 光栅器（`renderer.hpp`）改为两阶段：
  1. 三角形建立阶段按三角形并行；
  2. 光栅阶段按每 8 行一个行带并行。
  - 消除了原实现中重叠三角形之间的数据竞争，结果与线程数无关、逐位一致。
  - 包围盒为 NaN 的三角形直接剔除。
- 着色函数改为 `const`；`MooShader`、`MooMaterial`、`FBO`、`texture`、`Point_Light` 的成员都给了确定的初值。

### 健壮性（Legacy 头文件）
- `array/character/database/motion/nnet/robot/texture` 等文件中，`assert` 和不检查返回值的 `fread/fwrite`，改为 `MOO_CHECK_ALWAYS`：截断或损坏的二进制、打不开文件、内存不足时，所有构建配置都会给出明确报错。
- `mesh.hpp`：变长数组（VLA）改为 `std::vector`，MSVC 可以编译。
- `model.hpp`：对 malloc 出来的 `MooShader` 先就地构造再使用。
- 去掉了 UTF-8 BOM；`__restrict__` 改为跨编译器的 `MOO_RESTRICT`。

### 运行与验证
- 帧率节流：优先用高精度可等待计时器，不可用时退回 1ms 精度的 sleep。
- 控制台输出 `Cost`、`Interval` 和 `Gpu`（GPU 时间戳）。
- 自动验证用的命令行参数：`--frames`、`--capture-frame`、`--renderer`、`--hidden`、`--unlimited-fps`。
- 新增工具 `MooImageDiff`（逐像素对拍，退出码分级）、`verify_gl.py`、`gl_probe.cpp`。

---

## 二、编辑器与运动控制（本轮）

### 1. 目录与模块重组
- 根目录的 25 个 `.hpp` 和 `main.cpp` 移到 `Source/`，按 UE 风格分为 Public/Private：

  | 模块 | 内容 |
  |---|---|
  | `Runtime/Core` | INTERFACE 库：`Core.h`、`GlContext.h`、`ViewportTexture.h`；不依赖 GL、窗口或旧数学头 |
  | `Runtime/RHI` | `GlRhi`：GL 对象 RAII、时间戳、着色器读取 |
  | `Runtime/Renderer` | `RenderTypes`、`GpuRenderer`、`CpuImagePresenter` |
  | `Runtime/ApplicationCore` | `Platform`（XInput、帧节流）、`AppWindow`（GLFW、ImGui）、`Input`（键鼠采集与路由） |
  | `Sandbox` | `SceneTypes`、`SandboxScene`（PImpl 门面），以及 Private 下的 `SceneRender.inl`、`MotionMatching.inl` |
  | `Editor` | `EditorUi`（停靠式编辑器） |
  | `Launch` | `Main.cpp`（6 行）、`Application`（窗口、帧节拍、输入路由、截图组装） |
  | `Legacy` | 原有的数学/动画/仿真/CPU 渲染头文件 |

- GLSL 移到 `Shaders/`，工具移到 `Tools/ImageDiff`、`Tools/Verify`。
- CMake 拆成 `MooCore`、`MooRHI`、`MooRenderer`、`MooApplicationCore`、`MooSandbox`、`MooEditor`，再加 `MooLab.exe`。
  - 每个模块只导出 Public 包含路径。
  - 配置时检查公共头文件名在全仓唯一。
- 依赖方向：
  - Renderer → RHI/Core；ApplicationCore → Core；Sandbox → Renderer/Core。
  - Editor 只使用 Sandbox 的公共纯数据和 ImGui。
- Legacy 头里有大量非 inline 定义，所以只允许 `SandboxScene.cpp` 包含；两个 `.inl` 在 CMake 中标记为 `HEADER_FILE_ONLY`，避免 ODR 重复定义。

### 2. 场景封装与 main 瘦身
- 原先 main 里的场景逻辑收进 `FSandboxScene`（PImpl）。对外只暴露 `SceneTypes.h` 中的纯数据：
  - `FSceneSettings`：可调参数；
  - `FSceneInput`：移动/视角输入；
  - `FSceneCommand`：模式切换、暂停、单步、重置、恢复默认；
  - `FSceneSnapshot`：只读状态，含 22 个骨骼视图和 MM 统计；
  - `FFrameStats`。
- 两种场景模式：
  - `Physics`：原刚体加力矩控制的演示，图像与历史基线逐字节一致。
  - `MotionMatching`：角色控制。
- 支持暂停、单步、重置和恢复默认。设置值非法时拒绝，并在 UI 中显示原因。
- 资源路径改为编译期宏，不再修改进程工作目录：
  - 资源目录不存在时，警告并退回到工作目录下的 `../resources/`；
  - 着色器打不开或编译失败时，明确退出。

### 3. 停靠式编辑器（ImGui Docking）
- 窗口 1600×900，布局如下：
  - 左侧：场景大纲（Camera、PointLight、Character 和 22 个胶囊）。
  - 中间：Viewport，离屏 720×720 图像等比留边显示；纹理无效时不绘制。
  - 右侧：属性面板。
  - 底部：播放/仿真、Motion Matching、统计三个页签。
- 支持"恢复布局"，不写 `imgui.ini`。中文标签使用系统已装字体，加载失败时退回英文。
- 可调参数：
  - 物理：重力、阻尼、子步、接触平面与 eps、力矩模式、控制截止步；
  - 相机：跟随距离、高度、灵敏度；
  - 灯光：是否跟随、位置、辐射基值/振幅；
  - 显示：胶囊 roughness，胶囊/接触点/角色皮肤开关；
  - MM：速度上限、匹配权重。
- Metallic 控件置灰，并注明旧着色模型不使用它（本次不改 BRDF）。

### 4. 键鼠输入
- 新增 `FInputRouter`：
  - WASD 按相机的水平朝向控制角色（S 是转身向后走，不是倒退）。
  - 方向键环绕相机；空格暂停/继续；`.` 单步。
- 键盘只在窗口和 Viewport 都聚焦、且没有输入框占用时生效；鼠标只操作 UI。
- 键盘和 XInput 手柄的输入先相加，再做径向死区。
- 失焦时整体清空输入。门控关闭期间也记录按键状态，避免在控件里按住的键切回视口后被当成新按下。
- 窗口开启粘滞键，修复了快捷键短按丢失的问题。

### 5. Motion Matching 重写（`MotionMatching.inl`）
- 只读 `resources/long_motion.bin` 和 `database.bin`，在内存中搜索和混合，不改写资源，也不调用会影响 Physics 基线的旧拼接/对齐函数。
- 候选筛选（共 53500 帧）：
  - 搜索窗口 `[s-1, s+75]` 不能跨越片段接缝。
  - 含倒走帧的窗口剔除，因此 S 键只会转身行走。
  - 存储的速度和脚部特征要与重算值一致（容差分别为 0.1 m/s、0.02 m），且必须是有限值、不超过速度和角速度上限。
  - 角速度特征在内存中重算，不使用数据库里的 avels。
- 搜索：全量确定性 argmin。代价 = 速度、偏航角速度、脚部位置三项归一化后加权求和；按接触状态过滤；代价相同时取最小源帧。
- 触发：每 50 步强制搜索一次，并带 10 步去抖。
- 播放与混合：
  - 最多 5 条播放轨；新轨用 25 步 smoothstep 淡入，完全接管后丢弃更旧的轨。
  - 根旋转拆成偏航和残差两部分；偏航用不回绕的连续值做标量插值，避免轨间偏航差超过 180° 时绕最短弧翻转。
- 期望速度和偏航用临界阻尼弹簧平滑（半衰期 0.15 s / 0.2 s）；相机前向退化时退回 +z。
- 数据校验失败时禁止进入 MM，Physics 仍可正常使用。

### 6. 验证工具与文档
- 新增 `Tools/Verify/SceneProbe.cpp`，检查以下内容：
  - Physics 重置、模式往返后的历史帧；
  - 暂停与单步；
  - MM 边界情况；
  - 3 万步长测与运动预算。
  - 退出码：0 通过，1 硬失败，3 预算失败。
- `verify_gl.py` 新增 CLI 负例。
- 命令行新增 `--scene physics|motion-matching`。
- README 同步了构建、目录、编辑器、键盘、命令行参数、验证脚本和性能数据。
- `.qoder/agents` 中 6 个专家定义同步了新目录和资源路径的约定。

### 7. 评审修复（爱莎两轮，最终 0 阻断）
- 资源加载统一走 `ResourceRoot()`，不再切换工作目录。
- 构造时检查胶囊数与骨骼数一致（包括 MatchShapes）；快照循环以 `Bones.size()` 为上限。
- `ApplyCommand`：
  - 拒绝未知 Mode；
  - "恢复默认"加了检查；
  - 改为"设置模式或重置时才执行 Reset"。
- 视角输入过滤非有限值。
- MM 的 `Reset` 和 `Tick` 增加前置检查。
- 大纲条目数改为由快照决定。
- `GlRhi.cpp` 在宏未定义时回退到 `Shaders/`。
- 补充显式的 `#include <tuple>`，并修正几处注释中的旧文件名。

---

## 三、验证结果

> 本节记录的是脚部 IK 落地**之前**的状态。`foot_contacts.bin` 重新烘焙后角色姿态已变，下表的截图哈希与性能数字已被第五节取代；当前基准见 README（CPU 60/300 帧 `F32E92295DFD7FFF` / `6CC141F39195243E`）。

| 项目 | 结果 |
|---|---|
| 构建 | Release、Debug 均成功（`MooLab`、`SceneProbe`、`MooImageDiff`、`gl_probe`） |
| 截图哈希 | 与历史一致（已被第五节取代）：GPU 60 帧 `2490D9D5…2BE4`，GPU 300 帧 `E451D89D…C94D`；CPU 60/300 帧 `004F883938E6AD0A` / `EB3E20DCBC32FF56` |
| SceneProbe | Release、Debug 均退出码 0 |
| `verify_gl.py probe` | 通过 |
| 编码与 diff-check | 通过 |
| 手动 UI 测试 | 用户确认通过（方向键环绕相机、WASD 长按） |

**Motion Matching 运动预算**（Release）：

| 场景 | 结果 |
|---|---|
| 按 W 起步 | 第 24 步到位，方向误差 16.8° |
| 松键停步 | 第 50 步停下 |
| 90° 转向 | 第 44 步完成，误差 3.4° |
| 180° 转身 | 第 47 步完成，转过 176.7°；混合期间有 13 步倒走 |
| 过渡跳变 | 单步位移最大 0.072 m，偏航最大 7.37° |
| 3 万步长测 | 结果确定可复现；每次搜索约 0.162 ms，每步约 0.177 ms |

**性能**（Release，60 帧节流，3×600 帧）：
- GPU 路径 `Cost` P50 约 2.91 ms，在历史基线 2.83–2.94 ms 范围内。
- CPU 路径约 8.40 ms，历史基线为 8.73–8.99 ms。

---

## 四、待定事项
- 180° 转身混合期间的 13 步倒走，二选一：
  - (a) 调整验收口径；
  - (b) 改为在局部速度空间混合。
- `Source/Legacy` 中 14 个文件有内容改动（见第一部分），提交时是单独拆一个提交，还是合并后在说明里列出。
- 构建环境方面的遗留项：
  - `VCPKG_ROOT` 已写入 User 级环境变量（2.0 提交前完成）；
  - `compile_gl_probe.cmd` 里写死了 VS 2026 的 vcpkg/CMake 路径；
  - `.qoder/`（已迁移至仓库根目录）与 `run_ui_case.py` 已随 2.0 入库；
  - 基线放在被 git 忽略的目录里，且 Pillow 依赖没有声明。
- 仓库路径须为 ASCII（爱莎登记项）。
- `.qoder/repowiki` 原基于 `fadced7` 生成，已按 2.0 现状手工更新（`_index.yaml`、依赖/构建/路径卡片重写，`zh/content` 深度文档加版本横幅并重写三篇顶层文档）；后续仍建议在 IDE 中重新生成全量知识库。

---

## 五、逐脚触地标注与 MM 脚锁 IK（本轮，未验收）

### 改动
- 新增 `MooAnimation` 静态库：
  - `FootContacts`：逐脚触地标注，以及 `foot_contacts.bin` 的读写和校验；
  - `TwoBoneIk`：稳定两骨 IK、脚尖支点抬跟、soft 伸展。
- 新增资源 `resources/foot_contacts.bin`（53564 B，SHA256 `CB345AE3…E0B3`），绑定 `long_motion.bin` 的 FNV-1a64 `78313e50ce582fa2`。旧资源未改。
- 新增工具 `FootContactBake`：用 `BakeFootContacts` 重新生成标签。Release、Debug 两版生成的结果与资源逐字节相同。
- Motion Matching：
  - 检索改为按 2 bit 逐脚标签精确过滤，删除了原来的脚尖速度支撑推断；
  - 标签缺失或损坏时 MM 不可用，Physics 照常运行。
- `FootIk.inl`：显示姿态的脚锁后处理。
  - 只改两腿的 UpLeg/Leg/Foot 旋转，以及膝、踝、趾的位置；
  - 脚趾保持动画的世界旋转；
  - 不回馈 Matcher，root 和骨盆不变。
- 已批准的参数：抬跟最多 45°、soft 区 5 cm、膝角每步最多 14.5°，几何冲突时膝角优先。释放条件为关闭、标签下降、水平距离超过 0.3 m、残差超过 2 cm。
- 设置 `bFootIk` 默认开启；编辑器的 MM 面板显示每只脚的标签、锁定状态和残差；快照新增 `Feet`。
- 用户批准的淡出风险修复：
  - 权重归零后的交还判据改为 float 世界坐标下膝、踝、趾的米制位置误差，容差为 `max(1e-5, 32*FLT_EPSILON*坐标尺度)`，且不得有几何冲突；匹配时直接交还本步原始动画，不保留重解残差。
  - 增加 30 步终止兜底，以及 `MaxSettleSteps` / `ForcedHandbacks` 快照诊断；探针要求强制交还次数为 0，并监测最长连续淡出。
  - 增加真实行走至原点外 100 m 后关闭 IK 的回归检查与 one-ULP 脚尖偏移负对照。
  - UI 的 `clamp` 标记改为实际残差，并标出膝角优先冲突。
- 与计划的偏差：
  - 标签最短触地 5 帧、最长空隙 2 帧，脚尖使用水平速度；
  - 锁点只固定水平位置；
  - 膝弯参考轴从源动作标定；
  - 新增抬跟，Foot 不再保持动画旋转；
  - 抬跟另留 20° 数值余量。
- `verify_gl.py`：证据目录默认改为 `build/vs2026/foot_ik_v1/`，`ground_pbr_v1/` 只读；`input_hashes` 纳入新标签和烘焙工具。

### 验证结果
| 项目 | 结果 |
|---|---|
| 构建 | 淡出修复后 Debug、Release 六个目标均成功（MooLab、MooImageDiff、SceneProbe、FootContactBake、gl_probe、MaterialProbe）；均只有存量 C4244 警告 |
| 标注 | 左/右触地率 0.545/0.545；镜像一致率 ≥ 0.99992；静止帧双支撑 0.938；9 类损坏文件全部被拒绝 |
| 硬检查 | 淡出修复后 Release 默认线程、Release 单线程、Debug 完整 SceneProbe 均为 `hard checks=PASS`，`budgetFailures=4`（退出码 3）；两版独立 `--foot` 结果相同 |
| 淡出 | 行走 2013 步至 100.036 m 后关闭 IK，11 步内交还；固定最长连续淡出 24/24 步，随机 37/39 步；最长 settle 1/1 步，强制交还 0/0 次 |
| 运动摘要 | `physics600=a0690f9545f3dd1c`（不变）；修复后 `matching3000=3c276a0a51a8d80e`，三种配置及地面编辑后一致；修复前为 `b222084b59d6705f` |
| 历史证据 | 303 个文件 SHA256 不变；Physics 7 配置截图及三配置 Physics/hidden 截图均与 ground_pbr_v1 逐字节相同 |
| MM 截图 | 三配置 60/300 帧 CPU/GPU 按原阈值对拍通过；全部 8 类场景截图跨配置逐字节一致 |
| 脚锁预算 | 膝步 14.50°、滑移比 2.9% / 18.5%、锁定覆盖率通过；**4 项失败**（见下） |
| 原 MM 预算 | 单步位移最大 0.074 m，偏航最大 9.77°，均通过 |
| 性能（修复前独立测量，Release） | 标签加载加哈希约 50 ms；IK 每步平均约 0.001–0.002 ms，最大 0.19 ms；本次验证含并行任务，计时不用于性能优劣结论 |

4 项预算失败（阈值没有放宽）：
- 固定脚本：锁定期间单步滑移 2/692 步超过 1 mm，最大 13.9 mm。
- 随机脚本：锁定期间单步滑移 37/26271 步超过 1 mm，最大 14.1 mm。
- 固定脚本：释放后附加位移 8 次超过 20 mm，最大 54.0 mm。
- 随机脚本：释放后附加位移 141 次超过 20 mm，最大 66.3 mm。

### 证据与评审
- 本次日志位于 `build/vs2026/foot_ik_v1/probe_tmp/fade_*`；最新 Debug 完整结果为 `fade_scene_debug_retry.log`。第一次完整 Debug 在 900 秒看门狗处超时，重跑仅将工具等待时间增加到 1800 秒，没有改变验收阈值，最终退出 3。
- 两版 `FootContactBake` 再次写入隔离副本，SHA256 均为 `CB345AE34E7FF895EBA8359DBE578707788250281347717DA843D0632303E0B3`，与资源相同。
- 已重新执行 `capture`：7 配置、材质、gl_probe 通过；随后被 `scene_release` 退出 3 正确拦下。只完成部分候选验证，未生成合格场景收据，未冻结基线。
- 用户授权的追加编码评审：本次五文件增量 0 红旗、3 建议；第 2 轮的淡出不结束风险已修补，并有运行证据。

### 待定
- 4 项预算仍失败，数值未改善，因此**没有冻结新基线**，`build/vs2026/foot_ik_v1/baseline/` 不存在；整体功能仍未验收。
- 追加评审的覆盖建议尚未全部实施：30 步强制交还边界未受控触发；100 m 关闭用例实际覆盖右脚锁定，没有分别从左右脚远点锁定状态关闭，也没有全腿动画对照的共用负对照谓词。
- 数值容差随世界坐标增长；100 m 结果不证明任意远距离的误差都不可见，30 步兜底也不是移动动画目标的收敛证明。
- 手动 UI 验收未进行。
