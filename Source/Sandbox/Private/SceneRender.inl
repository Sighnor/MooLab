// 仅由SandboxScene.cpp包含；保留默认绘制顺序和旧接触点投影，地面先于全部模型绘制。
void FSandboxScene::InitializeGpu(FGpuRenderer& Renderer)
{
    MOO_CHECK_ALWAYS(Impl->GpuItems.empty(), "scene GPU resources already initialized");
    std::map<const vec3*, FTextureId> Cache;
    auto Upload = [&](const texture& Image)
    {
        const auto Found = Cache.find(Image.colors);
        if(Found != Cache.end()) return Found->second;
        const FTextureId Id = Renderer.CreateTextureRGBF(ImageView(Image));
        Cache.emplace(Image.colors, Id);
        return Id;
    };
    for(size_t Index = 0; Index < Impl->CpuRenderer.models.size(); ++Index)
    {
        MooModel* Model = Impl->CpuRenderer.models[Index];
        for(int MeshIndex = 0; MeshIndex < Model->mesh_count; ++MeshIndex)
        {
            FGpuDrawItem Item;
            Item.Model = Model;
            Item.ModelIndex = static_cast<int>(Index);
            Item.Mesh = Renderer.CreateMesh(MeshView(Model->meshes[MeshIndex]));
            const int ShaderType = Model->shaders[MeshIndex].shader_type;
            MOO_CHECK_ALWAYS(ShaderType == LIGHT || ShaderType == PBR, "unsupported GPU shader type %d", ShaderType);
            if(ShaderType == PBR)
            {
                Item.Shading = EShading::Pbr;
                const MooMaterial& Material = Model->materials[MeshIndex];
                Item.Pbr.Roughness = Material.roughness;
                Item.Pbr.Metallic = Material.metallic;
                Item.Pbr.AlbedoMap = Upload(Material.tex);
            }
            Impl->GpuItems.push_back(Item);
        }
    }
    // 材质标量始终以Settings为准，首次GPU渲染即反映当前配置
    Impl->ApplyMaterials();
}
FFrameParams FSandboxScene::FImpl::MakeFrameParams()
{
    FFrameParams Frame;
    CopyRenderMatrix(Frame.View, Camera.get_view_matrix());
    CopyRenderMatrix(Frame.Proj, Camera.get_projection_matrix());
    std::memcpy(Frame.CameraPos, &Camera.position, sizeof(vec3));
    std::memcpy(Frame.LightPos, &Light.light_position, sizeof(vec3));
    std::memcpy(Frame.LightRadiance, &Light.light_radiance, sizeof(vec3));
    Frame.ClearColor[0] = 100.f / 255.f;
    Frame.ClearColor[1] = 100.f / 255.f;
    Frame.ClearColor[2] = 200.f / 255.f;
    Frame.ZNear = Camera.z_near;
    Frame.ClearDepth = 100.f;
    return Frame;
}
bool FSandboxScene::FImpl::MakeGroundParams(const FFrameParams& Frame, int Width, int Height, FGroundParams& Out) const
{
    if(!Settings.Ground.bVisible) return false;
    const std::array<double, 4> Plane = GetGroundPlane();
    FGroundParams Result;
    Result.Roughness = Settings.Ground.Material.Roughness;
    Result.Metallic = Settings.Ground.Material.Metallic;
    if(!BuildGroundParams(Frame, Plane.data(), Settings.Ground.CellSize, Width, Height, Result)) return false;
    Out = Result;
    return true;
}
void FSandboxScene::FImpl::DrawGroundCpu(const FFrameParams& FrameParams, FBO& Target)
{
    FGroundParams Ground;
    if(!MakeGroundParams(FrameParams, Target.cols, Target.rows, Ground)) return;
    const vec3 Normal(Ground.Plane[0], Ground.Plane[1], Ground.Plane[2]);
    const vec3 CameraPos = Camera.position;
    const vec3 LightPos = Light.light_position;
    const vec3 LightRadiance = Light.light_radiance;
    // 按行并行，每个像素只写自身颜色与深度
#pragma omp parallel for schedule(static)
    for(int Row = 0; Row < Target.rows; ++Row)
    {
        for(int Column = 0; Column < Target.cols; ++Column)
        {
            FGroundSample Sample;
            if(!EvaluateGroundPixel(Ground, FrameParams.CameraPos, FrameParams.ZNear, FrameParams.ClearDepth, Column, Row, Sample)) continue;
            const vec3 Albedo(Sample.Albedo[0], Sample.Albedo[1], Sample.Albedo[2]);
            const vec3 Position(Sample.Position[0], Sample.Position[1], Sample.Position[2]);
            const vec3 Color = MooShadePbr(Albedo, Normal, Position, CameraPos, LightPos, LightRadiance, Ground.Roughness, Ground.Metallic);
            Target.getcolor(Row, Column) = clampv(255.f * Color);
            Target.getdepth(Row, Column) = Sample.Depth;
        }
    }
}
void FSandboxScene::FImpl::DrawContacts(FGpuRenderer* Renderer)
{
    if(!Settings.bShowContacts) return;
    std::vector<int> Pixels;
    for(int Index = 1; Index < Motion.nbones(); ++Index)
    {
        if(!Shapes().bone_shapes(Index - 1).ifcontact) continue;
        if(Renderer)
        {
            int X = 0, Y = 0;
            if(ProjectPointToPixel(Positions(Index), Camera, Camera.fbo[0].rows, Camera.fbo[0].cols, X, Y))
            {
                Pixels.push_back(X);
                Pixels.push_back(Y);
            }
        }
        else draw_point(Positions(Index), Camera, &Camera.fbo[0], vec3(0.f, 255.f, 0.f));
    }
    if(Renderer)
    {
        const float Color[3] = {0.f, 1.f, 0.f};
        Renderer->DrawPixelSquares(Pixels.data(), static_cast<int>(Pixels.size() / 2), 2, Color);
    }
}
void FSandboxScene::BuildRenderData(FGpuRenderer& Renderer)
{
    const FFrameParams Frame = Impl->MakeFrameParams();
    if(Impl->Settings.bShowCharacter)
    {
        const auto View = MeshView(Impl->SkinMesh);
        for(const auto& Item : Impl->GpuItems)
            if(Item.ModelIndex == CharacterModelIndex) Renderer.UpdateMeshVertices(Item.Mesh, View.Positions, View.Normals, View.VertexCount);
    }
    Renderer.BeginFrame(Frame);
    // 地面先于模型，写入的view depth参与模型深度测试；与CPU同一尺寸（fbo）准备参数
    FGroundParams Ground;
    if(Impl->MakeGroundParams(Frame, Impl->Camera.fbo[0].cols, Impl->Camera.fbo[0].rows, Ground)) Renderer.DrawGround(Ground);
    for(const auto& Item : Impl->GpuItems)
    {
        if(Item.ModelIndex == CharacterModelIndex && !Impl->Settings.bShowCharacter) continue;
        if(Item.ModelIndex >= FirstCapsuleModelIndex && !Impl->Settings.bShowCapsules) continue;
        float Matrix[16];
        CopyRenderMatrix(Matrix, Item.Model->transform);
        Renderer.DrawMesh(Item.Mesh, Matrix, Item.Shading, Item.Shading == EShading::Pbr ? &Item.Pbr : nullptr);
    }
    Impl->DrawContacts(&Renderer);
    Renderer.EndFrame();
}
FImageRGBFView FSandboxScene::RenderCpu()
{
    auto& State = *Impl;
    auto& Frame = State.Camera.fbo[0];
    Frame.set(vec3(100.f, 100.f, 200.f), 100.f);
    State.DrawGroundCpu(State.MakeFrameParams(), Frame);
    draw(&State.LightModel, State.Camera, &Frame, &State.Parameters, true);
    if(State.Settings.bShowCharacter) draw(&State.SkinModel, State.Camera, &Frame, &State.Parameters, true);
    if(State.Settings.bShowCapsules)
        for(int Index = 0; Index < State.CapsuleModels.size; ++Index)
            draw(&State.CapsuleModels(Index), State.Camera, &Frame, &State.Parameters, true);
    State.DrawContacts(nullptr);
    return {reinterpret_cast<const float*>(Frame.colors), Frame.cols, Frame.rows};
}
void FSandboxScene::ReadCpuRGB8(FImageRGB8& Image)
{
    // 使用旧OpenCV饱和舍入，输出通用RGB；应用层独立执行RGB->BGR写图。
    auto& Frame = Impl->Camera.fbo[0];
    cv::Mat FloatImage(Frame.rows, Frame.cols, CV_32FC3, Frame.colors), Bytes;
    FloatImage.convertTo(Bytes, CV_8UC3, 1.0f);
    Image.Width = Frame.cols;
    Image.Height = Frame.rows;
    Image.Pixels.assign(Bytes.data, Bytes.data + Bytes.total() * Bytes.elemSize());
}
