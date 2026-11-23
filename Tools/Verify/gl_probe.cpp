#include "AppWindow.h"
#include "GpuRenderer.h"
#include "PbrShading.h"
#include "CpuImagePresenter.h"
#include "GlRhi.h"
#include "Core.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

// 临时验收探针，仅包含新接口，不进入生产目标。
struct FInvalidContext final : IGlContextOwner
{
    bool IsGlContextCurrent() const override { return false; }
};

int main(int Argc, char** Argv)
{
#if defined(_MSC_VER)
    // 工具负例须无交互退出；只关闭本探针进程的CRT abort弹窗，不改生产错误处理
    _set_abort_behavior(0, _WRITE_ABORT_MSG | _CALL_REPORTFAULT);
#endif
    if(Argc > 1 && std::strcmp(Argv[1], "bad-context") == 0)
    {
        FInvalidContext Invalid;
        FGpuRenderer Renderer(Invalid, 31, 17);
        return 0;
    }
    FAppWindow Window({31, 17, "MOOLAB probe", true});
    FGpuRenderer Renderer(Window, 31, 17);
    FFrameParams Frame;
    float Model[16] = {};
    for(int Index = 0; Index < 16; Index += 5)
    {
        Model[Index] = Frame.View[Index] = Frame.Proj[Index] = 1.f;
    }
    Frame.ClearColor[0] = Frame.ClearColor[1] = 100.f / 255.f;
    Frame.ClearColor[2] = 200.f / 255.f;
    Frame.LightRadiance[0] = 1.f;
    Frame.ZNear = 0.05f;
    Frame.ClearDepth = 128.f;
    FImageRGB8 Image;
    Renderer.BeginFrame(Frame);
    Renderer.EndFrame();
    // 故意污染DSA仍会读取的像素传输状态，Readback必须自行清理。
    FGlBuffer PixelBuffer = CreateGlBuffer(65536, nullptr, GL_DYNAMIC_STORAGE_BIT);
    glBindBuffer(GL_PIXEL_PACK_BUFFER, PixelBuffer.Get());
    glPixelStorei(GL_PACK_ROW_LENGTH, 43);
    glPixelStorei(GL_PACK_SKIP_ROWS, 2);
    glPixelStorei(GL_PACK_SKIP_PIXELS, 3);
    Renderer.ReadbackRGB8(Image);
    MOO_CHECK_ALWAYS(Image.Width == 31 && Image.Height == 17, "rectangular target size wrong");
    for(size_t Pixel = 0; Pixel < Image.Pixels.size(); Pixel += 3)
    {
        MOO_CHECK_ALWAYS(Image.Pixels[Pixel] == 100 && Image.Pixels[Pixel + 1] == 100 && Image.Pixels[Pixel + 2] == 200, "empty scene clear wrong");
    }
    std::puts("empty rectangular / polluted pack: PASS");

    float Positions[] = {-0.75f,-0.75f,-0.5f, 0.75f,-0.75f,-0.5f, 0.f,0.75f,-0.5f};
    const float Normals[] = {0,0,1, 0,0,1, 0,0,1};
    const float Uvs[] = {0,0, 0,1, 1,1};
    const uint16_t Indices[] = {0,1,2};
    if(Argc > 1 && std::strcmp(Argv[1], "bad-index") == 0)
    {
        const uint16_t BadIndices[] = {0,1,3};
        Renderer.CreateMesh({Positions, Normals, Uvs, BadIndices, 3, 1});
        return 0;
    }
    if(Argc > 1 && std::strcmp(Argv[1], "bad-present") == 0)
    {
        FCpuImagePresenter Presenter(Window, 31, 17);
        std::vector<float> Colors(30 * 17 * 3, 0.f);
        Presenter.Present({Colors.data(), 30, 17}, 31, 17);
        return 0;
    }
    FMeshId Mesh = Renderer.CreateMesh({Positions, Normals, Uvs, Indices, 3, 1});
    if(Argc > 1 && std::strcmp(Argv[1], "bad-update") == 0)
    {
        Renderer.UpdateMeshVertices(Mesh, Positions, Normals, 2);
        return 0;
    }
    auto Draw = [&]()
    {
        Renderer.BeginFrame(Frame);
        Renderer.DrawMesh(Mesh, Model, EShading::Light, nullptr);
        Renderer.EndFrame();
        Renderer.ReadbackRGB8(Image);
    };
    auto CountRed = [&]()
    {
        int Count = 0;
        for(size_t Pixel = 0; Pixel < Image.Pixels.size(); Pixel += 3)
        {
            Count += Image.Pixels[Pixel] == 255 && Image.Pixels[Pixel + 1] == 0 && Image.Pixels[Pixel + 2] == 0;
        }
        return Count;
    };
    Draw();
    MOO_CHECK_ALWAYS(CountRed() > 0, "light triangle invisible");
    const auto Before = Image.Pixels;
    Renderer.UpdateMeshVertices(Mesh, Positions, Normals, 3);
    Draw();
    MOO_CHECK_ALWAYS(Before == Image.Pixels, "identity vertex update changed image");
    for(int Index = 0; Index < 3; Index++) { Positions[Index * 3] += 0.25f; }
    Renderer.UpdateMeshVertices(Mesh, Positions, Normals, 3);
    Draw();
    MOO_CHECK_ALWAYS(Before != Image.Pixels, "changed vertices did not affect image");
    for(float Depth : {0.05f, 128.f, 129.f})
    {
        for(int Index = 0; Index < 3; Index++) { Positions[Index * 3 + 2] = -Depth; }
        Renderer.UpdateMeshVertices(Mesh, Positions, Normals, 3);
        Draw();
        MOO_CHECK_ALWAYS((CountRed() > 0) == (Depth == 128.f), "near/clear depth comparison failed at %f", Depth);
    }
    std::puts("vertex updates / near boundary / power-of-two clear depth: PASS");

    // 程序化地面：相机在原点、单位View/Proj，n=(0,0,1)、d=2 的正对平面对全部射线 t=2（恒深）
    auto MakeGround = [&](const double (&Plane)[4])
    {
        FGroundParams Ground;
        MOO_CHECK_ALWAYS(BuildGroundParams(Frame, Plane, 0.5f, 31, 17, Ground), "ground params rejected");
        return Ground;
    };
    Frame.ClearDepth = 100.f;
    const double Facing[4] = {0., 0., 1., 2.};
    const FGroundParams FacingGround = MakeGround(Facing);
    if(Argc > 1 && std::strcmp(Argv[1], "bad-ground-outside") == 0)
    {
        Renderer.DrawGround(FacingGround);
        return 0;
    }
    if(Argc > 1 && std::strncmp(Argv[1], "bad-ground-", 11) == 0)
    {
        Renderer.BeginFrame(Frame);
        if(std::strcmp(Argv[1], "bad-ground-size") == 0)
        {
            FGroundParams WrongSize = FacingGround;
            WrongSize.TargetHeight++;
            Renderer.DrawGround(WrongSize);
            return 0;
        }
        if(std::strcmp(Argv[1], "bad-ground-after-mesh") == 0) Renderer.DrawMesh(Mesh, Model, EShading::Light, nullptr);
        if(std::strcmp(Argv[1], "bad-ground-after-points") == 0)
        {
            const int Point[] = {5, 5};
            const float Color[] = {0, 1, 0};
            Renderer.DrawPixelSquares(Point, 1, 2, Color);
        }
        if(std::strcmp(Argv[1], "bad-ground-twice") == 0) Renderer.DrawGround(FacingGround);
        Renderer.DrawGround(FacingGround);
        return 0;
    }
    auto IsClear = [&](int Row, int Col)
    {
        const size_t Pixel = static_cast<size_t>(Row * 31 + Col) * 3;
        return Image.Pixels[Pixel] == 100 && Image.Pixels[Pixel + 1] == 100 && Image.Pixels[Pixel + 2] == 200;
    };
    auto DrawGround = [&](const FGroundParams& Ground, bool bMesh)
    {
        Renderer.BeginFrame(Frame);
        Renderer.DrawGround(Ground);
        if(bMesh) Renderer.DrawMesh(Mesh, Model, EShading::Light, nullptr);
        Renderer.EndFrame();
        Renderer.ReadbackRGB8(Image);
    };
    // GPU覆盖掩码必须与CPU求交逐像素一致，返回覆盖像素数
    auto CheckMask = [&](const FGroundParams& Ground, const char* Name)
    {
        DrawGround(Ground, false);
        int Covered = 0;
        for(int Row = 0; Row < 17; Row++)
        {
            for(int Col = 0; Col < 31; Col++)
            {
                FGroundSample Sample;
                const bool bCpu = EvaluateGroundPixel(Ground, Frame.CameraPos, Frame.ZNear, Frame.ClearDepth, Col, Row, Sample);
                MOO_CHECK_ALWAYS(bCpu == !IsClear(Row, Col), "%s ground mask differs from CPU at row %d col %d", Name, Row, Col);
                Covered += bCpu;
            }
        }
        return Covered;
    };
    MOO_CHECK_ALWAYS(CheckMask(FacingGround, "facing") == 31 * 17, "facing ground must cover last column and bottom row");
    const double Horizon[4] = {0., 1., 0., 0.5};
    MOO_CHECK_ALWAYS(CheckMask(MakeGround(Horizon), "horizon") == 9 * 31, "horizontal ground must cover rows 8..16 only");
    const double Slanted[4] = {1., 1., 0., 0.5};
    const int SlantedCovered = CheckMask(MakeGround(Slanted), "slanted");
    MOO_CHECK_ALWAYS(SlantedCovered > 0 && SlantedCovered < 31 * 17, "slanted ground must partially cover");
    const double AtNear[4] = {0., 0., 1., 0.05};
    const double PastNear[4] = {0., 0., 1., 0.0501};
    const double AtClear[4] = {0., 0., 1., 100.};
    const double PastClear[4] = {0., 0., 1., 100.01};
    MOO_CHECK_ALWAYS(CheckMask(MakeGround(AtNear), "near") == 0 && CheckMask(MakeGround(PastNear), "past-near") == 31 * 17,
                     "ground t == ZNear rejected, beyond accepted");
    MOO_CHECK_ALWAYS(CheckMask(MakeGround(AtClear), "clear") == 31 * 17 && CheckMask(MakeGround(PastClear), "past-clear") == 0,
                     "ground t == ClearDepth accepted, beyond rejected");
    // 绕过BuildGroundParams的相机侧检查，着色器自身也须拒绝背面
    FGroundParams Backface = FacingGround;
    for(float& Value : Backface.Plane) Value = -Value;
    MOO_CHECK_ALWAYS(CheckMask(Backface, "backface") == 0, "camera behind ground drawn");
    std::puts("ground mask vs CPU / edge row+column / near+clear / horizon / slant / backface: PASS");

    // 地面写入view depth=2：更近的恒深三角形可见，更远的被遮挡且画面与仅地面一致
    DrawGround(FacingGround, false);
    const auto GroundOnly = Image.Pixels;
    for(float Depth : {1.9f, 2.1f})
    {
        const float Triangle[] = {-0.75f,-0.75f,-Depth, 0.75f,-0.75f,-Depth, 0.f,0.75f,-Depth};
        Renderer.UpdateMeshVertices(Mesh, Triangle, Normals, 3);
        DrawGround(FacingGround, true);
        MOO_CHECK_ALWAYS((CountRed() > 0) == (Depth < 2.f), "ground depth occlusion failed at %f", Depth);
        MOO_CHECK_ALWAYS((Image.Pixels == GroundOnly) == (Depth > 2.f), "occluded mesh changed ground image at %f", Depth);
    }
    std::puts("ground depth occlusion with constant-depth mesh: PASS");

    Renderer.BeginFrame(Frame);
    const int Points[] = {0,0,30,16};
    const float Green[] = {0,1,0};
    Renderer.DrawPixelSquares(Points, 2, 2, Green);
    Renderer.EndFrame();
    Renderer.ReadbackRGB8(Image);
    int GreenCount = 0;
    for(int Row = 0; Row < 17; Row++)
    {
        for(int Col = 0; Col < 31; Col++)
        {
            const size_t Pixel = static_cast<size_t>(Row * 31 + Col) * 3;
            const bool bGreen = Image.Pixels[Pixel] == 0 && Image.Pixels[Pixel + 1] == 255 && Image.Pixels[Pixel + 2] == 0;
            MOO_CHECK_ALWAYS(bGreen == ((Row < 3 && Col < 3) || (Row >= 14 && Col >= 28)), "edge point mask wrong");
            GreenCount += bGreen;
        }
    }
    MOO_CHECK_ALWAYS(GreenCount == 18, "edge square clipping wrong");
    std::puts("edge pixel squares: PASS");

    FCpuImagePresenter Presenter(Window, 31, 17);
    std::vector<float> Colors(31 * 17 * 3, 0.f);
    for(size_t Pixel = 0; Pixel < Colors.size(); Pixel += 3) { Colors[Pixel] = 255.f; }
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, PixelBuffer.Get());
    glPixelStorei(GL_UNPACK_ROW_LENGTH, 43);
    glPixelStorei(GL_UNPACK_SKIP_ROWS, 2);
    glPixelStorei(GL_UNPACK_SKIP_PIXELS, 3);
    glPixelStorei(GL_UNPACK_SWAP_BYTES, GL_TRUE);
    Presenter.Present({Colors.data(), 31, 17}, 31, 17);
    std::vector<uint8_t> CpuImage(31 * 17 * 3);
    glReadBuffer(GL_BACK);
    glReadPixels(0, 0, 31, 17, GL_RGB, GL_UNSIGNED_BYTE, CpuImage.data());
    for(size_t Pixel = 0; Pixel < CpuImage.size(); Pixel += 3)
    {
        MOO_CHECK_ALWAYS(CpuImage[Pixel] == 255 && CpuImage[Pixel + 1] == 0 && CpuImage[Pixel + 2] == 0, "polluted unpack presenter wrong");
    }
    std::puts("polluted unpack / CPU presentation: PASS");
    FGlGpuTimer Timer;
    for(int Index = 0; Index < 2000; Index++) { Timer.Begin(); Timer.End(); }
    glFlush();
    MOO_CHECK_ALWAYS(glGetError() == GL_NO_ERROR, "GL error in probes");
    std::puts("timer ring stress / GL error check: PASS");
    return 0;
}
