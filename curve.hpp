#ifndef MOOLAB_CURVE
#define MOOLAB_CURVE

#include "array.hpp"
#include "global.hpp"
#include "vec.hpp"

template<typename T>
T recursive_bezier(const slice1d<T> points, float t) 
{
    array1d<T> sub_points(points.size - 1);

    if(points.size == 2)
    {
        return (1 - t) * points(0) + t * points(1);
    }
    
    for(int i = 0; i < (points.size - 1); i++)
    {
       sub_points(i) = ((1 - t) * points(i) + t * points(i + 1));
    }

    return recursive_bezier(sub_points,t);
}

template<typename T>
array1d<T> bezier_curve(const slice1d<T> points, int sample_num = 1000) 
{
    array1d<T> curve_points(sample_num + 1);

    for(float t = 0.f; t <= 1.f; t += 1.f / sample_num) 
    {
        T point = recursive_bezier(points, t);
        curve_points(t * sample_num) = point;
    }

    return curve_points;
}

template<typename T>
void pid_control(
        T &x,
        T &vel,
        T &acc,
        T x_goal,
        T p_const, 
        T i_const, 
        T d_const, 
        float k, 
        float dt)
{
    static T ek_0 = T(0.f);
    static T ek_1 = T(0.f);
    static T ek_2 = T(0.f);

    ek_2 = ek_1;
    ek_1 = ek_0;
    ek_0 = x_goal - x;

    acc = k * ((p_const * (ek_0 - ek_1) + i_const * ek_0 - d_const * (ek_0 - 2 * ek_1 + ek_2)));
    vel = vel + acc * dt;
    x = x + vel * dt;
}

template<typename T>
array1d<T> PID_curve(const slice1d<T> points, int sample_num = 300)
{
    array1d<T> curve_points(sample_num + 1);

    T x = points(0);
    T vel = T(0.f);
    T acc = T(0.f);

    for(float t = 0; t < 1.f; t += 1.f / sample_num)
    {
        float i = t * (points.size - 1);
        int i0 = floor(i);
        int i1 = ceil(i);
        T point = lerp(points(i0), points(i1), i - i0);

        pid_control(
            x,
            vel,
            acc,
            point,
            T(1.f), 
            T(0.05f), 
            T(0.2f), 
            200.f, 
            1.f / 60.f);

        curve_points(t * sample_num) = x;
    }

    return curve_points;
}

template<typename T>
void cubic_spine_interpolation(
        T &x,
        const T &x0,
        const T &x1,
        const T &v0,
        const T &v1,
        float t)
{
    T a = 2 * x0 - 2 * x1 + v0 + v1;
    T b = - 3 * x0 + 3 * x1 - 2 * v0 - v1;
    T c = v0;
    T d = x0;

    x = a * t * t * t + b * t * t + c * t + d;
}

template<typename T>
array1d<T> spine_curve(const slice1d<T> points, float mSampleTime = 0.05f)
{
    array1d<T> curve_points((points.size - 1) / mSampleTime + 1);
    array1d<T> vec(points.size);

    vec(0) = points(1) - points(0);
    vec(points.size - 1) = points(points.size - 1) - points(points.size - 2);
    for(int i = 1; i < points.size - 1; i++)
    {
        vec(i) = (points(i + 1) - points(i - 1)) / 2.f;
    }

    int id = 0;

    for(float t = 0.f; t < 1.f; t += 1.f / (curve_points.size - 1), id++)
    {
        float i = t * (points.size - 1);
        int i0 = floor(i);
        int i1 = ceil(i);

        T point;

        cubic_spine_interpolation(
            point, 
            points(i0), 
            points(i1), 
            vec(i0), 
            vec(i0), 
            i - i0);

        curve_points(id) = point;
    }

    curve_points(curve_points.size - 1) = points(points.size - 1);

    return curve_points;
}

void LFPB_function0(
        float &y,
        float y0,
        float vel,
        float acc0,
        float t0,
        float tb0,
        float t)
{
    y = y0 + (vel + 0.5 * acc0 * t - acc0 * (t0 + tb0)) * t + (-vel + acc0 * (0.5 * t0 + tb0)) * t0;
}

void LFPB_function1(
        float &y,
        float y0,
        float vel,
        float acc0,
        float t0,
        float tb0,
        float t)
{
    y = y0 + vel * (t - t0) - 0.5 * acc0 * tb0 * tb0;
}

void LFPB_function2(
        float &y,
        float y1,
        float vel,
        float acc1,
        float t1,
        float tb1,
        float t)
{
    y = y1 + (vel + 0.5 * acc1 * t - acc1 * (t1 - tb1)) * t + (-vel + acc1 * (0.5 * t1 - tb1)) * t1;
}

array1d<float> LFPB(const slice1d<float> points, const slice1d<float> td, float acc_value, float mSampleTime = 0.05f)
{
    array1d<float> t(points.size);
    array1d<float> tb(points.size);
    array1d<float> acc(points.size);
    array1d<float> vel(points.size);

    t(0) = 0;
    for (int i = 1;i < points.size; i++)
    {
        t(i) = t(i - 1) + td(i - 1);
    }

    array1d<float> curve_points(t(points.size - 1) / mSampleTime + 1);

    //初始加速度和结束加速度的方向
    if (points(1) == points(0))
    {
        acc(0) = 0;
        tb(0) = 0;
    }
    else
    {
        acc(0) = sgn(points(1) - points(0)) * acc_value;
        tb(0) = td(0) - sqrt(td(0) * td(0) - 2 * (points(1) - points(0)) / acc(0));
    }
    if (points(points.size - 2) == points(points.size - 1))
    {
        acc(points.size - 1) = 0;
        tb(points.size - 1) = 0;
    }
    else
    {
        acc(points.size - 1) = sgn(points(points.size - 2) - points(points.size - 1)) * acc_value;
        tb(points.size - 1) = td(points.size - 2) - sqrt(td(points.size - 2) * td(points.size - 2) - 2 * (points(points.size - 2) - points(points.size - 1)) / acc(points.size - 1));
    }
    //第一段直线速度和最后一段直线速度
    vel(0) = (points(1) - points(0)) / (td(0) - 0.5 * tb(0));
    vel((points.size - 2)) = (points(points.size - 1) - points(points.size - 2)) / (td(points.size - 2) - 0.5 * tb(points.size - 1));
    //计算剩下每段的直线速度
    for (int i = 1;i < points.size - 2;i++)
    {
        vel(i) = (points(i + 1) - points(i)) / td(i);
    }
    //计算剩下每段的加速度和加速时间
    for (int i = 1;i < points.size - 1;i++)
    {
        if (vel(i) == vel((i - 1)))
        {
            acc(i) = 0;
            tb(i) = 0;
        }
        else
        {
            acc(i) = sgn(vel(i) - vel(i - 1)) * acc_value;
            tb(i) = (vel(i) - vel(i - 1)) / (2 * acc(i));
        }
    }
    // //位移补偿量
    for (int k = 1;k < points.size - 1;k++)
    {
        points(k) = points(k) + 0.5 * acc(k) * tb(k) * tb(k);
    }
    //计算每个时刻的角度
    for (int k = 0;k < points.size - 1;k++)
    {
        //ceil向上取整，保证大于对应函数范围下界
        for (int i = ceil(t(k) / mSampleTime);i < (t(k) + tb(k)) / mSampleTime;i++)
        {
            LFPB_function0(curve_points(i), points(k), vel(k), acc(k), t(k), tb(k), i * mSampleTime);
        }
        for (int i = ceil((t(k) + tb(k)) / mSampleTime);i < (t(k + 1) - tb((k + 1))) / mSampleTime;i++)
        {
            LFPB_function1(curve_points(i), points(k), vel(k), acc(k), t(k), tb(k), i * mSampleTime);
        }
        for (int i = ceil((t(k + 1) - tb((k + 1))) / mSampleTime);i < t(k + 1) / mSampleTime;i++)
        {
            LFPB_function2(curve_points(i), points(k + 1), vel(k), acc(k + 1), t(k + 1), tb(k + 1), i * mSampleTime);
        }
    }
   
    return curve_points;
}

#endif