#include "TwoBoneIk.h"

#include <cmath>
#include <limits>

namespace
{
using FVec = std::array<double, 3>;
using FQuat = std::array<double, 4>;
constexpr double Pi = 3.14159265358979323846;

FVec Add(const FVec& A, const FVec& B) { return {A[0] + B[0], A[1] + B[1], A[2] + B[2]}; }
FVec Sub(const FVec& A, const FVec& B) { return {A[0] - B[0], A[1] - B[1], A[2] - B[2]}; }
FVec Scale(const FVec& A, double S) { return {A[0] * S, A[1] * S, A[2] * S}; }
double Dot(const FVec& A, const FVec& B) { return A[0] * B[0] + A[1] * B[1] + A[2] * B[2]; }
FVec Cross(const FVec& A, const FVec& B)
{
    return {A[1] * B[2] - A[2] * B[1], A[2] * B[0] - A[0] * B[2], A[0] * B[1] - A[1] * B[0]};
}
double Length(const FVec& A) { return std::sqrt(Dot(A, A)); }
double Clamp(double Value, double Min, double Max) { return Value < Min ? Min : (Value > Max ? Max : Value); }
bool IsFinite(const FVec& A) { return std::isfinite(A[0]) && std::isfinite(A[1]) && std::isfinite(A[2]); }
bool IsFinite(const FQuat& Q) { return std::isfinite(Q[0]) && std::isfinite(Q[1]) && std::isfinite(Q[2]) && std::isfinite(Q[3]); }

FQuat Multiply(const FQuat& A, const FQuat& B)
{
    return {A[0] * B[0] - A[1] * B[1] - A[2] * B[2] - A[3] * B[3],
            A[0] * B[1] + A[1] * B[0] + A[2] * B[3] - A[3] * B[2],
            A[0] * B[2] - A[1] * B[3] + A[2] * B[0] + A[3] * B[1],
            A[0] * B[3] + A[1] * B[2] - A[2] * B[1] + A[3] * B[0]};
}
bool Normalize(FQuat& Q)
{
    const double Norm = std::sqrt(Q[0] * Q[0] + Q[1] * Q[1] + Q[2] * Q[2] + Q[3] * Q[3]);
    if(!(Norm > 0.) || !std::isfinite(Norm)) return false;
    for(double& Value : Q) Value /= Norm;
    return true;
}
FVec Rotate(const FQuat& Q, const FVec& V)
{
    // v + w*t + u x t，t = 2 u x v
    const FVec U = {Q[1], Q[2], Q[3]};
    const FVec T = Scale(Cross(U, V), 2.);
    return Add(Add(V, Scale(T, Q[0])), Cross(U, T));
}
// 绕单位轴Axis转Angle弧度
FQuat AxisAngle(const FVec& Axis, double Angle)
{
    const double Half = 0.5 * Angle;
    const double S = std::sin(Half);
    return {std::cos(Half), Axis[0] * S, Axis[1] * S, Axis[2] * S};
}
// 绕Axis从From到To的有向角
double SignedAngle(const FVec& From, const FVec& To, const FVec& Axis)
{
    return std::atan2(Dot(Axis, Cross(From, To)), Dot(From, To));
}
// 与Direction垂直的单位向量：Candidate去掉沿Direction分量；退化时改用与Direction最不平行的世界轴。
FVec Perpendicular(const FVec& Direction, const FVec& Candidate, bool& bFallback)
{
    FVec Result = Sub(Candidate, Scale(Direction, Dot(Candidate, Direction)));
    double Norm = Length(Result);
    bFallback = !(Norm > 1e-6 * (Length(Candidate) + 1e-300));
    if(bFallback)
    {
        const double Ax = std::fabs(Direction[0]), Ay = std::fabs(Direction[1]), Az = std::fabs(Direction[2]);
        const FVec Axis = Ax <= Ay && Ax <= Az ? FVec{1., 0., 0.} : (Ay <= Az ? FVec{0., 1., 0.} : FVec{0., 0., 1.});
        Result = Sub(Axis, Scale(Direction, Dot(Axis, Direction)));
        Norm = Length(Result);
    }
    return Scale(Result, 1. / Norm);
}
}

bool SolveTwoBoneIk(const FTwoBoneIkInput& Input, FTwoBoneIkOutput& Out)
{
    FQuat Upper = Input.UpperRotation;
    FQuat Lower = Input.LowerRotation;
    const FVec& A = Input.Root;
    Out = FTwoBoneIkOutput{};
    const double L1 = Length(Input.MidOffset);
    const double L2 = Length(Input.EndOffset);
    const bool bValid = IsFinite(A) && IsFinite(Upper) && IsFinite(Lower) && IsFinite(Input.MidOffset) && IsFinite(Input.EndOffset)
        && IsFinite(Input.Target) && IsFinite(Input.ReferenceAxis) && Normalize(Upper) && Normalize(Lower) && L1 >= 1e-6 && L2 >= 1e-6;
    if(!bValid)
    {
        Out.UpperRotation = Input.UpperRotation;
        Out.LowerRotation = Input.LowerRotation;
        Out.Mid = Add(A, Rotate(Input.UpperRotation, Input.MidOffset));
        Out.End = Add(Out.Mid, Rotate(Input.LowerRotation, Input.EndOffset));
        return false;
    }
    const FVec B = Add(A, Rotate(Upper, Input.MidOffset));
    const FVec C = Add(B, Rotate(Lower, Input.EndOffset));
    const FVec AB = Sub(B, A), AC = Sub(C, A), ToTarget = Sub(Input.Target, A);
    const double CurrentReach = Length(AC);
    const double RawDistance = Length(ToTarget);
    const double MinDistance = std::fabs(L1 - L2) + 1e-3 * (L1 + L2);
    const double MaxDistance = std::fmin(L1 + L2, std::fmax(0.999 * (L1 + L2), CurrentReach));
    const double Distance = Clamp(RawDistance, MinDistance, MaxDistance);
    Out.Unreachable = std::fabs(RawDistance - Distance);

    // 弯曲轴
    const FVec DirectionAB = Scale(AB, 1. / L1);
    const FVec Normal = Cross(AC, AB);
    const double NormalLength = Length(Normal);
    FVec Axis;
    if(NormalLength > 1e-3 * CurrentReach * L1)
    {
        Axis = Scale(Normal, (Dot(Normal, Input.ReferenceAxis) < 0. ? -1. : 1.) / NormalLength);
    }
    else
    {
        bool bFallback = false;
        Axis = Perpendicular(DirectionAB, Input.ReferenceAxis, bFallback);
        Out.bDegenerateAxis = true;
    }

    // 当前有向角与余弦定理目标角；直腿时a->c退化为a->b方向
    const FVec DirectionAC = CurrentReach > 1e-9 ? Scale(AC, 1. / CurrentReach) : DirectionAB;
    const double Alpha0 = SignedAngle(DirectionAC, DirectionAB, Axis);
    double Beta0 = SignedAngle(Sub(A, B), Sub(C, B), Axis);
    if(Beta0 < 0.) Beta0 += 2. * Pi;
    const double Alpha = std::acos(Clamp((L1 * L1 + Distance * Distance - L2 * L2) / (2. * L1 * Distance), -1., 1.));
    const double Beta = std::acos(Clamp((L1 * L1 + L2 * L2 - Distance * Distance) / (2. * L1 * L2), -1., 1.));
    const double Theta1 = Alpha - Alpha0;
    const double Theta2 = Beta - Beta0;
    Upper = Multiply(AxisAngle(Axis, Theta1), Upper);
    Lower = Multiply(AxisAngle(Axis, Theta1 + Theta2), Lower);
    Normalize(Upper);
    Normalize(Lower);

    // 整链摆向目标
    const FVec Reach = Add(Rotate(Upper, Input.MidOffset), Rotate(Lower, Input.EndOffset));
    const double ReachLength = Length(Reach);
    const FVec From = ReachLength > 1e-12 ? Scale(Reach, 1. / ReachLength) : DirectionAC;
    const FVec To = RawDistance > 1e-9 ? Scale(ToTarget, 1. / RawDistance) : DirectionAC;
    const double W = 1. + Dot(From, To);
    FQuat Swing;
    if(W < 1e-6)
    {
        bool bFallback = false;
        const FVec Flip = Perpendicular(From, Axis, bFallback);
        Swing = {0., Flip[0], Flip[1], Flip[2]};
    }
    else
    {
        const FVec V = Cross(From, To);
        Swing = {W, V[0], V[1], V[2]};
        Normalize(Swing);
    }
    Upper = Multiply(Swing, Upper);
    Lower = Multiply(Swing, Lower);
    Normalize(Upper);
    Normalize(Lower);

    Out.UpperRotation = Upper;
    Out.LowerRotation = Lower;
    Out.Mid = Add(A, Rotate(Upper, Input.MidOffset));
    Out.End = Add(Out.Mid, Rotate(Lower, Input.EndOffset));
    return true;
}

double SoftenFootReach(double D, double MaxDistance, double Zone)
{
    if(!std::isfinite(D) || !std::isfinite(MaxDistance) || !std::isfinite(Zone) || D < 0. || MaxDistance < 0. || Zone <= 0.)
        return std::numeric_limits<double>::quiet_NaN();
    Zone = std::fmin(Zone, MaxDistance);
    if(Zone == 0.) return 0.;
    const double Start = MaxDistance - Zone;
    return D <= Start ? D : MaxDistance - Zone * std::exp(-(D - Start) / Zone);
}

bool SolveFootPivotIk(const FFootPivotIkInput& Input, FFootPivotIkOutput& Out)
{
    constexpr double MaxLift = Pi / 4.;
    constexpr double LiftMargin = 20. * Pi / 180.;
    constexpr double KneeStep = 14.5 * Pi / 180.;
    Out = FFootPivotIkOutput{};
    FTwoBoneIkInput Chain = Input.Chain;
    FQuat Foot = Input.FootRotation;
    const double L1 = Length(Chain.MidOffset), L2 = Length(Chain.EndOffset), FootLength = Length(Input.ToeOffset);
    if(!IsFinite(Chain.Root) || !IsFinite(Chain.UpperRotation) || !IsFinite(Chain.LowerRotation) || !IsFinite(Chain.MidOffset)
        || !IsFinite(Chain.EndOffset) || !IsFinite(Chain.ReferenceAxis) || !IsFinite(Foot) || !IsFinite(Input.ToeOffset)
        || !IsFinite(Input.ToeTarget) || !IsFinite(Input.PreviousLiftAxis)
        || !std::isfinite(Input.PreviousKneeAngle) || Input.PreviousKneeAngle > Pi
        || !std::isfinite(Input.SoftWeight) || Input.SoftWeight < 0. || Input.SoftWeight > 1.
        || !Normalize(Chain.UpperRotation) || !Normalize(Chain.LowerRotation) || !Normalize(Foot)
        || !std::isfinite(L1) || !std::isfinite(L2) || !std::isfinite(FootLength)
        || L1 < 1e-6 || L2 < 1e-6 || FootLength < 1e-6) return false;

    const FVec AnimatedMid = Add(Chain.Root, Rotate(Chain.UpperRotation, Chain.MidOffset));
    const FVec AnimatedEnd = Add(AnimatedMid, Rotate(Chain.LowerRotation, Chain.EndOffset));
    const double AnimatedDistance = Length(Sub(AnimatedEnd, Chain.Root));
    const double MinDistance = std::fabs(L1 - L2) + 1e-3 * (L1 + L2);
    const double MaxDistance = std::fmin(L1 + L2, std::fmax(0.999 * (L1 + L2), AnimatedDistance));
    if(MinDistance >= 0.999 * (L1 + L2)) return false; // 极端骨长比不能满足稳定伸展预留量。
    auto KneeAt = [&](double D) { return std::acos(Clamp((L1 * L1 + L2 * L2 - D * D) / (2. * L1 * L2), -1., 1.)); };
    auto DistanceAtKnee = [&](double B) { return std::sqrt(std::fmax(0., L1 * L1 + L2 * L2 - 2. * L1 * L2 * std::cos(B))); };
    const double PreviousKnee = Input.PreviousKneeAngle < 0. ? KneeAt(AnimatedDistance) : Input.PreviousKneeAngle;
    const double BandMin = std::fmax(MinDistance, DistanceAtKnee(Clamp(PreviousKnee - KneeStep, 0., Pi)));
    const double BandMax = std::fmin(MaxDistance, DistanceAtKnee(Clamp(PreviousKnee + KneeStep, 0., Pi)));
    if(BandMin > BandMax + 1e-12) return false;

    // u是脚尖到踝的向量；从同一旋转与偏移构造目标及输出，避免混用float动画位置导致锁点残差。
    const FVec U = Scale(Rotate(Foot, Input.ToeOffset), -1.);
    const FVec W = Sub(Chain.Root, Input.ToeTarget);
    const FVec Forward = {-U[0], 0., -U[2]};
    FVec Axis = Cross({0., 1., 0.}, Forward);
    double AxisLength = Length(Axis);
    const bool bVertical = AxisLength <= 1e-3 * FootLength;
    if(bVertical)
    {
        Axis = {Input.PreviousLiftAxis[0], 0., Input.PreviousLiftAxis[2]};
        AxisLength = Length(Axis);
        if(AxisLength <= 1e-9)
        {
            Axis = Cross({0., 1., 0.}, {-W[0], 0., -W[2]});
            AxisLength = Length(Axis);
        }
        if(AxisLength <= 1e-9) { Axis = {1., 0., 0.}; AxisLength = 1.; }
    }
    Axis = Scale(Axis, 1. / AxisLength);
    if(!bVertical && Cross(Axis, U)[1] < 0.) Axis = Scale(Axis, -1.); // 竖直退化保留历史符号，不按数值噪声翻轴。
    Out.LiftAxis = Axis;
    const FVec UPerp = Sub(U, Scale(Axis, Dot(U, Axis)));
    const FVec WPerp = Sub(W, Scale(Axis, Dot(W, Axis)));
    const double Rho = Length(UPerp), Omega = Length(WPerp), Product = Rho * Omega;
    const double Along = Dot(Sub(U, W), Axis);
    const double Constant = Along * Along + Rho * Rho + Omega * Omega;
    const double Psi = std::atan2(Dot(Axis, Cross(UPerp, WPerp)), Dot(UPerp, WPerp));
    constexpr double Base = 0.;
    const double Limit = Product > 1e-10 ? std::fmax(Base, std::fmin(MaxLift, Psi - LiftMargin)) : Base;
    auto ReachAt = [&](double Angle) { return Length(Sub(Add(Input.ToeTarget, Rotate(AxisAngle(Axis, Angle), U)), Chain.Root)); };
    const double BaseDistance = ReachAt(Base), LiftDistance = ReachAt(Limit);
    const double SoftMax = 0.999 * (L1 + L2);
    const double Soft = SoftenFootReach(BaseDistance, SoftMax, std::fmin(0.05, 0.5 * (SoftMax - MinDistance)));
    const double Reference = std::fmin(BaseDistance, std::fmax(AnimatedDistance, Soft));
    const double Desired = BaseDistance - Input.SoftWeight * (BaseDistance - Reference);
    const double Low = std::fmax(BandMin, LiftDistance), High = std::fmin(BandMax, BaseDistance);
    Out.bConflict = Low > High + 1e-12;
    double ActualDistance;
    if(!Out.bConflict) ActualDistance = Clamp(Desired, Low, std::fmax(Low, High));
    else ActualDistance = LiftDistance > BandMax ? BandMax : BandMin; // 取膝角带中离几何区间最近的端点。

    double Lift = Base;
    if(ActualDistance < LiftDistance) Lift = Limit;
    else if(ActualDistance < BaseDistance - 1e-12 && Product > 1e-10 && Limit > Base)
    {
        const double Cosine = Clamp((Constant - ActualDistance * ActualDistance) / (2. * Product), -1., 1.);
        Lift = Clamp(Psi - std::acos(Cosine), Base, Limit);
    }
    const FQuat LiftRotation = AxisAngle(Axis, Lift);
    const FVec AnkleRequest = Add(Input.ToeTarget, Rotate(LiftRotation, U));
    const FVec ToAnkle = Sub(AnkleRequest, Chain.Root);
    const double RequestedDistance = Length(ToAnkle);
    const FVec Direction = RequestedDistance > 1e-9 ? Scale(ToAnkle, 1. / RequestedDistance)
        : (AnimatedDistance > 1e-9 ? Scale(Sub(AnimatedEnd, Chain.Root), 1. / AnimatedDistance)
                                 : Scale(Rotate(Chain.UpperRotation, Chain.MidOffset), 1. / L1));
    // 无冲突时直接使用几何目标；冲突时沿该方向落在膝角带内，残差由最终FK测得。
    Chain.Target = Out.bConflict ? Add(Chain.Root, Scale(Direction, ActualDistance)) : AnkleRequest;
    if(!SolveTwoBoneIk(Chain, Out.Chain)) return false;
    Out.FootRotation = Multiply(LiftRotation, Foot);
    Normalize(Out.FootRotation);
    Out.Toe = Add(Out.Chain.End, Rotate(Out.FootRotation, Input.ToeOffset));
    Out.LiftAngle = Lift;
    Out.KneeAngle = KneeAt(Length(Sub(Out.Chain.End, Chain.Root)));
    Out.Residual = Length(Sub(Out.Toe, Input.ToeTarget));
    return IsFinite(Out.Toe) && IsFinite(Out.FootRotation) && std::isfinite(Out.Residual) && std::isfinite(Out.KneeAngle);
}
