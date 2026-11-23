#pragma once

// 上下文契约不暴露GL或窗口类型；GL资源必须早于上下文析构。
class IGlContextOwner
{
public:
    virtual bool IsGlContextCurrent() const = 0;
protected:
    ~IGlContextOwner() = default;
};
