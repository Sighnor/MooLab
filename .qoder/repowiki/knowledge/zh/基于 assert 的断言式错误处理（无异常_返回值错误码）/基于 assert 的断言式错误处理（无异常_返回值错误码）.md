---
kind: error_handling
name: 基于 assert 的断言式错误处理（无异常/返回值错误码）
category: error_handling
scope:
    - '**'
source_files:
    - array.hpp
    - character.hpp
    - database.hpp
    - motion.hpp
    - nnet.hpp
    - main.cpp
---

## 1. 使用的系统/方法

该仓库是一个基于 raylib 的 C++ 演示程序，**没有使用任何结构化错误处理机制**：
- 未使用 `std::exception`、自定义异常类型或 `throw`。
- 未定义错误码枚举或错误结构体。
- 未使用 `try/catch` 块进行错误捕获。
- 未调用 `fprintf(stderr)` / `perror` / `exit` 等标准错误输出与终止 API。
- 唯一的错误检测手段是 C 标准库的 `assert`（`<cassert>`），用于在运行时检查前置条件。

## 2. 关键文件

- `array.hpp` — 最密集地使用 `assert` 的文件，覆盖一维/二维数组的边界检查、空指针检查、维度一致性校验。
- `character.hpp`、`database.hpp`、`motion.hpp`、`nnet.hpp` — 在文件 I/O 打开后使用 `assert(f != NULL)` 检查文件句柄。
- `main.cpp` — 仅有一处被注释掉的 `assert(fft != NULL);`，当前运行路径中不产生实际错误分支。

## 3. 架构与约定

### 约定：所有“错误”被视为不可恢复的内部不变量违反

代码将资源缺失、索引越界、维度不匹配等视为 **内部编程错误**，而非可恢复的外部错误。处理方式统一为 `assert(...)`，即：
- 开发/调试构建下触发 abort，附带文件名与行号；
- 发布构建下 `assert` 被编译掉，静默继续执行（这是 C/C++ `assert` 的标准语义）。

### 具体模式

1. **数组访问边界检查**（`array.hpp`）：
   - 单维 `operator()(i)`：`assert(i >= 0 && i < size)`。
   - 二维 `operator()(i, j)`：`assert(i >= 0 && i < rows && j >= 0 && j < cols)`。
   - 切片 `slice1d` 构造：`assert(i >= 0 && i < rows)`。

2. **空指针检查**（`array.hpp` 及多个 I/O 头文件）：
   - `assert(data != NULL)` 在 `array2d_read` / `array2d_write` 等函数入口。
   - `assert(f != NULL)` 在 `Database_load/save`、`Motion_load/save`、`NNet_load/save` 等文件操作前。

3. **维度/类型一致性检查**（`array.hpp`、`nnet.hpp`）：
   - `assert((int)num == arr.size)`、`assert(arr.cols == res.size)` 等。
   - `assert(input.size == input_layer.size && output.size == output_layer.size)` 等网络层尺寸校验。

4. **无返回值错误传播**：
   - 所有 I/O 与计算函数返回 `void`，不通过返回值传递错误状态；调用方依赖 `assert` 保证参数合法。

## 4. 约定与约束

- **描述性约定**：本仓库将所有外部输入非法、资源缺失、索引越界视为内部不变量违反，通过 `assert` 表达，不使用异常、错误码或日志输出。
- **约束来源**：C 标准库 `<cassert>` 的语义由编译器/运行时强制——调试模式下失败即 abort，发布模式下断言被移除。因此该策略对发布构建不提供任何保护，属于纯开发期防护。
- **观察到的范围**：除 `main.cpp` 中被注释的一行外，全仓未发现其他形式的错误处理（无 `throw`、无 `try/catch`、无错误枚举、无 stderr 输出）。