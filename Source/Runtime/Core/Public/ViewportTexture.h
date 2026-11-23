#pragma once
#include <cstdint>

// 非拥有的显示视图。有效期由渲染器管理；UI不得释放Handle。
struct FViewportTexture
{
    uint64_t Handle = 0;
    int Width = 0;
    int Height = 0;
    float U0 = 0.f;
    float V0 = 1.f;
    float U1 = 1.f;
    float V1 = 0.f;
};
