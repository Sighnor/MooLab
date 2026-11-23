#pragma once
#include <array>

// 稳定两骨IK（如 髋a -> 膝b -> 踝c）。纯数值接口：内部双精度，不依赖存量头或GL，无共享可写状态。
// 四元数顺序 (w, x, y, z)，旋转均为世界系全局旋转；b = a + Upper*MidOffset，c = b + Lower*EndOffset。
//
//  1. 距离钳制：D = clamp(|t-a|, |l1-l2| + 1e-3(l1+l2), min(l1+l2, max(0.999(l1+l2), |c-a|)))
//  2. 弯曲轴：nc = (c-a)x(b-a)；|nc| > 1e-3|c-a||b-a| 时 n = sign(nc·Ref)*nc/|nc|，
//     否则（直腿/退化）n = Ref 去掉沿 a->b 分量后归一化，Ref 也退化时取与 a->b 最不平行的世界轴构造。
//  3. 按余弦定理得内角 alpha(髋)、beta(膝)，以绕 n 的有向角 alpha0、beta0(取[0,2pi)) 为起点：
//     Upper' = R(n, alpha-alpha0)*Upper，Lower' = R(n, alpha-alpha0+beta-beta0)*Lower，a->c 方向保持不变。
//  4. 最短弧 R2 把 a->c' 转到 a->t（|t-a|过小时保持原方向；反向时绕 n 转 180 度），作用于两骨。
// 骨长严格保持；末端骨（足）的全局旋转由调用方保持不变。
struct FTwoBoneIkInput
{
    std::array<double, 3> Root = {};
    std::array<double, 4> UpperRotation = {1., 0., 0., 0.};
    std::array<double, 4> LowerRotation = {1., 0., 0., 0.};
    std::array<double, 3> MidOffset = {};
    std::array<double, 3> EndOffset = {};
    std::array<double, 3> Target = {};
    // 世界系弯曲参考轴：决定膝弯方向并作为直腿退化时的回退轴，无需单位长
    std::array<double, 3> ReferenceAxis = {1., 0., 0.};
};

struct FTwoBoneIkOutput
{
    std::array<double, 4> UpperRotation = {1., 0., 0., 0.};
    std::array<double, 4> LowerRotation = {1., 0., 0., 0.};
    std::array<double, 3> Mid = {};
    std::array<double, 3> End = {};
    // 原始目标距离与钳后距离之差（m），即目标不可达量
    double Unreachable = 0.;
    bool bDegenerateAxis = false;
};

// 输入非有限或骨长过短（< 1e-6 m）时返回false，Out为按输入旋转得到的原姿态。
bool SolveTwoBoneIk(const FTwoBoneIkInput& Input, FTwoBoneIkOutput& Out);

// 指数soft区只决定期望伸展距离；由绕脚尖抬跟实现，不直接拖动踝目标。
// D、MaxDistance须有限且>=0，Zone须有限且>0；非法输入返回NaN。
double SoftenFootReach(double D, double MaxDistance, double Zone = 0.05);

// 脚尖支点IK：仅向上抬跟（最多45度），soft区5cm，膝角每步最多14.5度。
// 几何距离区间与膝角连续区间不相交时膝角优先；残差如实输出，不伪造脚尖锁定。
// 为避开距离反解的极值病态，抬跟还保留距最佳靠近角20度的余量，实际上限可小于45度。
// 与两骨求解器采用相同伸展上下限，保留动画直腿不动点；足旋转可改，趾旋转由调用方保持动画值。
struct FFootPivotIkInput
{
    FTwoBoneIkInput Chain;
    std::array<double, 4> FootRotation = {1., 0., 0., 0.};
    std::array<double, 3> ToeOffset = {};
    std::array<double, 3> ToeTarget = {};
    // 脚掌近乎竖直时沿用的世界系抬跟轴；进一步退化时取确定性世界轴。
    std::array<double, 3> PreviousLiftAxis = {1., 0., 0.};
    // 弧度；<0表示无历史，使用当前动画膝角。
    double PreviousKneeAngle = -1.;
    // 淡出时的soft权重[0,1]，几何抬跟范围始终完整保留。
    double SoftWeight = 1.;
};
struct FFootPivotIkOutput
{
    FTwoBoneIkOutput Chain;
    std::array<double, 4> FootRotation = {1., 0., 0., 0.};
    std::array<double, 3> Toe = {};
    std::array<double, 3> LiftAxis = {1., 0., 0.};
    double LiftAngle = 0.;
    double KneeAngle = 0.;
    // 期望脚尖与实际脚尖的距离，包含膝角优先产生的残差（m）。
    double Residual = 0.;
    bool bConflict = false;
};
// 非有限输入、非法状态、任一骨长过短或无法形成稳定距离区间时返回false；不含跨帧共享状态。
bool SolveFootPivotIk(const FFootPivotIkInput& Input, FFootPivotIkOutput& Out);
