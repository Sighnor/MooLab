// MooImageDiff：截图对拍工具（CPU参考 vs GPU输出），只依赖OpenCV imgcodecs。
//
// 用法：MooImageDiff A.png B.png [选项]
//   --shift-y N       以A(x, y)对比B(x, y - N)，只统计两图重叠的行，默认0（GPU路径已复刻CPU行映射，正常对拍不需要偏移）
//   --tol N           通道差（三通道最大值）<= N 视为“接近”，默认3
//   --gross N         通道差 > N 视为“严重偏差”，默认32
//   --min-close R     接近像素占比下限，默认0.99
//   --max-gross R     严重偏差像素占比上限，默认0.005
//   --diff-out PATH   输出差异图（每像素三通道最大差 x 8，灰度）
// 退出码：0 通过；1 未达阈值；2 参数或读图错误；3 尺寸不一致

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>

namespace
{
enum class EExitCode : int
{
    Pass = 0,
    ThresholdFailed = 1,
    UsageOrIoError = 2,
    SizeMismatch = 3,
};

int ToExitCode(EExitCode Code)
{
    return static_cast<int>(Code);
}

struct FDiffOptions
{
    std::string PathA;
    std::string PathB;
    std::string DiffOutPath;
    int ShiftY = 0;
    int Tolerance = 3;
    int GrossThreshold = 32;
    double MinCloseRatio = 0.99;
    double MaxGrossRatio = 0.005;
};

void PrintUsage()
{
    std::fprintf(stderr,
        "usage: MooImageDiff A.png B.png [--shift-y N] [--tol N] [--gross N] [--min-close R] [--max-gross R] [--diff-out PATH]\n"
        "exit: 0 pass, 1 threshold failed, 2 usage/io error, 3 size mismatch\n");
}

bool ParseInt(const char* Text, int& OutValue)
{
    char* End = nullptr;
    const long Value = std::strtol(Text, &End, 10);
    if(End == Text || *End != '\0' || Value < -100000 || Value > 100000)
    {
        return false;
    }
    OutValue = static_cast<int>(Value);
    return true;
}

bool ParseRatio(const char* Text, double& OutValue)
{
    char* End = nullptr;
    const double Value = std::strtod(Text, &End);
    if(End == Text || *End != '\0' || !(Value >= 0.0 && Value <= 1.0))
    {
        return false;
    }
    OutValue = Value;
    return true;
}

bool ParseOptions(int argc, char** argv, FDiffOptions& Options)
{
    int PositionalCount = 0;
    for(int ArgIndex = 1; ArgIndex < argc; ArgIndex++)
    {
        const char* Arg = argv[ArgIndex];
        const bool bHasValue = ArgIndex + 1 < argc;
        bool bOk = true;
        if(std::strcmp(Arg, "--shift-y") == 0 && bHasValue)
        {
            bOk = ParseInt(argv[++ArgIndex], Options.ShiftY);
        }
        else if(std::strcmp(Arg, "--tol") == 0 && bHasValue)
        {
            bOk = ParseInt(argv[++ArgIndex], Options.Tolerance) && Options.Tolerance >= 0;
        }
        else if(std::strcmp(Arg, "--gross") == 0 && bHasValue)
        {
            bOk = ParseInt(argv[++ArgIndex], Options.GrossThreshold) && Options.GrossThreshold >= 0;
        }
        else if(std::strcmp(Arg, "--min-close") == 0 && bHasValue)
        {
            bOk = ParseRatio(argv[++ArgIndex], Options.MinCloseRatio);
        }
        else if(std::strcmp(Arg, "--max-gross") == 0 && bHasValue)
        {
            bOk = ParseRatio(argv[++ArgIndex], Options.MaxGrossRatio);
        }
        else if(std::strcmp(Arg, "--diff-out") == 0 && bHasValue)
        {
            Options.DiffOutPath = argv[++ArgIndex];
        }
        else if(std::strncmp(Arg, "--", 2) != 0 && PositionalCount < 2)
        {
            (PositionalCount == 0 ? Options.PathA : Options.PathB) = Arg;
            PositionalCount++;
        }
        else
        {
            std::fprintf(stderr, "unknown or incomplete option: %s\n", Arg);
            return false;
        }
        if(!bOk)
        {
            std::fprintf(stderr, "invalid value for %s: %s\n", Arg, argv[ArgIndex]);
            return false;
        }
    }
    return PositionalCount == 2;
}
}

int main(int argc, char** argv)
{
    FDiffOptions Options;
    if(!ParseOptions(argc, argv, Options))
    {
        PrintUsage();
        return ToExitCode(EExitCode::UsageOrIoError);
    }

    const cv::Mat ImageA = cv::imread(Options.PathA, cv::IMREAD_COLOR);
    const cv::Mat ImageB = cv::imread(Options.PathB, cv::IMREAD_COLOR);
    if(ImageA.empty() || ImageB.empty())
    {
        std::fprintf(stderr, "failed to read image: %s\n", ImageA.empty() ? Options.PathA.c_str() : Options.PathB.c_str());
        return ToExitCode(EExitCode::UsageOrIoError);
    }
    if(ImageA.size() != ImageB.size())
    {
        std::fprintf(stderr, "size mismatch: A %dx%d, B %dx%d\n", ImageA.cols, ImageA.rows, ImageB.cols, ImageB.rows);
        return ToExitCode(EExitCode::SizeMismatch);
    }

    const int Width = ImageA.cols;
    const int Height = ImageA.rows;
    //A的第Y行对比B的第(Y - ShiftY)行，两者都须在图像范围内
    const int RowBegin = std::max(0, Options.ShiftY);
    const int RowEnd = std::min(Height, Height + Options.ShiftY);
    if(RowBegin >= RowEnd)
    {
        std::fprintf(stderr, "shift-y %d leaves no overlapping rows (height %d)\n", Options.ShiftY, Height);
        return ToExitCode(EExitCode::UsageOrIoError);
    }

    cv::Mat DiffImage;
    if(!Options.DiffOutPath.empty())
    {
        DiffImage = cv::Mat::zeros(Height, Width, CV_8UC1);
    }

    long long PixelCount = 0;
    long long IdenticalCount = 0;
    long long CloseCount = 0;
    long long GrossCount = 0;
    int MaxDiff = 0;
    double SquaredErrorSum = 0.0;
    for(int Y = RowBegin; Y < RowEnd; Y++)
    {
        const cv::Vec3b* RowA = ImageA.ptr<cv::Vec3b>(Y);
        const cv::Vec3b* RowB = ImageB.ptr<cv::Vec3b>(Y - Options.ShiftY);
        for(int X = 0; X < Width; X++)
        {
            int PixelDiff = 0;
            for(int Channel = 0; Channel < 3; Channel++)
            {
                const int ChannelDiff = std::abs(static_cast<int>(RowA[X][Channel]) - static_cast<int>(RowB[X][Channel]));
                PixelDiff = std::max(PixelDiff, ChannelDiff);
                SquaredErrorSum += static_cast<double>(ChannelDiff) * ChannelDiff;
            }
            PixelCount++;
            IdenticalCount += PixelDiff == 0 ? 1 : 0;
            CloseCount += PixelDiff <= Options.Tolerance ? 1 : 0;
            GrossCount += PixelDiff > Options.GrossThreshold ? 1 : 0;
            MaxDiff = std::max(MaxDiff, PixelDiff);
            if(!DiffImage.empty())
            {
                DiffImage.at<uint8_t>(Y, X) = static_cast<uint8_t>(std::min(255, PixelDiff * 8));
            }
        }
    }

    const double Mse = SquaredErrorSum / (static_cast<double>(PixelCount) * 3.0);
    const double CloseRatio = static_cast<double>(CloseCount) / PixelCount;
    const double GrossRatio = static_cast<double>(GrossCount) / PixelCount;
    const double IdenticalRatio = static_cast<double>(IdenticalCount) / PixelCount;
    const bool bPass = CloseRatio >= Options.MinCloseRatio && GrossRatio <= Options.MaxGrossRatio;

    char PsnrText[32];
    if(Mse > 0.0)
    {
        std::snprintf(PsnrText, sizeof(PsnrText), "%.2f", 10.0 * std::log10(255.0 * 255.0 / Mse));
    }
    else
    {
        std::snprintf(PsnrText, sizeof(PsnrText), "inf");
    }

    std::printf("size %dx%d shift-y %d rows %d..%d | max %d | psnr %s dB | identical %.4f%% | close(<=%d) %.4f%% (min %.4f%%) | gross(>%d) %.4f%% (max %.4f%%) | %s\n",
                Width, Height, Options.ShiftY, RowBegin, RowEnd - 1, MaxDiff, PsnrText,
                100.0 * IdenticalRatio,
                Options.Tolerance, 100.0 * CloseRatio, 100.0 * Options.MinCloseRatio,
                Options.GrossThreshold, 100.0 * GrossRatio, 100.0 * Options.MaxGrossRatio,
                bPass ? "PASS" : "FAIL");

    if(!DiffImage.empty() && !cv::imwrite(Options.DiffOutPath, DiffImage))
    {
        std::fprintf(stderr, "failed to write diff image: %s\n", Options.DiffOutPath.c_str());
        return ToExitCode(EExitCode::UsageOrIoError);
    }
    return ToExitCode(bPass ? EExitCode::Pass : EExitCode::ThresholdFailed);
}
