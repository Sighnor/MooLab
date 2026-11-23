#include "CpuImagePresenter.h"

#include "Core.h"
#include "GlRhi.h"
#include "RenderTypes.h"

struct FCpuImagePresenter::FImpl
{
    int Width = 0;
    int Height = 0;
    FGlProgram Program;
    FGlTexture2D Image;
    FGlTexture2D ViewportColor;
    FGlFramebuffer ViewportFramebuffer;
    //核心模式绘制必须绑定VAO；全屏三角形不读顶点属性，使用空VAO
    FGlVertexArray EmptyVertexArray;
};

FCpuImagePresenter::FCpuImagePresenter(const IGlContextOwner& Context, int Width, int Height)
{
    CheckGlContextReady(Context);
    Impl = std::make_unique<FImpl>();
    MOO_CHECK_ALWAYS(Width > 0 && Height > 0, "invalid image size %dx%d", Width, Height);
    Impl->Width = Width;
    Impl->Height = Height;
    Impl->Program = CreateGlProgramFromFiles(GetGlShaderPath("fullscreen.vert"), GetGlShaderPath("present.frag"));
    Impl->Image = CreateGlTexture2D(GL_RGB32F, Width, Height);
    Impl->ViewportColor = CreateGlTexture2D(GL_RGBA8, Width, Height);
    Impl->ViewportFramebuffer = CreateGlFramebuffer();
    glNamedFramebufferTexture(Impl->ViewportFramebuffer.Get(), GL_COLOR_ATTACHMENT0, Impl->ViewportColor.Get(), 0);
    MOO_CHECK_ALWAYS(glCheckNamedFramebufferStatus(Impl->ViewportFramebuffer.Get(), GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE,
                     "CPU viewport framebuffer incomplete");
    Impl->EmptyVertexArray = CreateGlVertexArray();
}

FCpuImagePresenter::~FCpuImagePresenter() = default;

void FCpuImagePresenter::UpdateViewport(const FImageRGBFView& Image)
{
    MOO_CHECK_ALWAYS(Image.Data != nullptr && Image.Width == Impl->Width && Image.Height == Impl->Height,
                     "CPU image does not match presenter");
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
    glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
    glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
    glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
    glPixelStorei(GL_UNPACK_SWAP_BYTES, GL_FALSE);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
    glTextureSubImage2D(Impl->Image.Get(), 0, 0, 0, Impl->Width, Impl->Height, GL_RGB, GL_FLOAT, Image.Data);
    glBindFramebuffer(GL_FRAMEBUFFER, Impl->ViewportFramebuffer.Get());
    glViewport(0, 0, Impl->Width, Impl->Height);
    glDisable(GL_SCISSOR_TEST);
    glDisable(GL_BLEND);
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_FRAMEBUFFER_SRGB);
    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
    glUseProgram(Impl->Program.Get());
    glBindTextureUnit(0, Impl->Image.Get());
    glBindVertexArray(Impl->EmptyVertexArray.Get());
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}
FViewportTexture FCpuImagePresenter::GetViewportTexture() const
{
    return {Impl->ViewportColor.Get(), Impl->Width, Impl->Height, 0.f, 1.f, 1.f, 0.f};
}

void FCpuImagePresenter::Present(const FImageRGBFView& Image, int BackbufferWidth, int BackbufferHeight)
{
    MOO_CHECK_ALWAYS(Image.Data != nullptr && Image.Width == Impl->Width && Image.Height == Impl->Height,
                     "CPU image %dx%d does not match presenter %dx%d", Image.Width, Image.Height, Impl->Width, Impl->Height);
    //最小化时帧缓冲尺寸为0，跳过呈现
    if(BackbufferWidth <= 0 || BackbufferHeight <= 0)
    {
        return;
    }

    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
    glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
    glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
    glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
    glPixelStorei(GL_UNPACK_SWAP_BYTES, GL_FALSE);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
    glTextureSubImage2D(Impl->Image.Get(), 0, 0, 0, Impl->Width, Impl->Height, GL_RGB, GL_FLOAT, Image.Data);

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glViewport(0, 0, BackbufferWidth, BackbufferHeight);
    glDisable(GL_SCISSOR_TEST);
    glDisable(GL_BLEND);
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);

    glUseProgram(Impl->Program.Get());
    glBindTextureUnit(0, Impl->Image.Get());
    glBindVertexArray(Impl->EmptyVertexArray.Get());
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}
