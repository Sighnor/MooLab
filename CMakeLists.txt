cmake_minimum_required(VERSION 3.21)
project(MooLab LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(OpenCV CONFIG REQUIRED COMPONENTS core imgproc imgcodecs)
find_package(OpenMP REQUIRED COMPONENTS CXX)
find_package(glfw3 CONFIG REQUIRED)
find_package(glad CONFIG REQUIRED)
find_package(imgui CONFIG REQUIRED)

add_library(MooBuildOptions INTERFACE)
if(MSVC)
    target_compile_definitions(MooBuildOptions INTERFACE NOMINMAX WIN32_LEAN_AND_MEAN _CRT_SECURE_NO_WARNINGS)
    target_compile_options(MooBuildOptions INTERFACE /utf-8 /permissive- /W3 /MP)
endif()
add_library(MooCore INTERFACE)
target_include_directories(MooCore INTERFACE Source/Runtime/Core/Public)
target_link_libraries(MooCore INTERFACE MooBuildOptions)

function(moo_module Target Directory)
    add_library(${Target} STATIC ${ARGN})
    target_include_directories(${Target} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/${Directory}/Public")
    target_link_libraries(${Target} PUBLIC MooCore)
    source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}/Source" FILES ${ARGN})
endfunction()

moo_module(MooRHI Source/Runtime/RHI
    Source/Runtime/RHI/Public/GlRhi.h Source/Runtime/RHI/Private/GlRhi.cpp)
target_link_libraries(MooRHI PRIVATE glad::glad)
target_compile_definitions(MooRHI PRIVATE MOOLAB_SHADER_DIR="${CMAKE_CURRENT_SOURCE_DIR}/Shaders")

moo_module(MooRenderer Source/Runtime/Renderer
    Source/Runtime/Renderer/Public/RenderTypes.h
    Source/Runtime/Renderer/Public/PbrShading.h Source/Runtime/Renderer/Private/PbrShading.cpp
    Source/Runtime/Renderer/Public/GpuRenderer.h Source/Runtime/Renderer/Private/GpuRenderer.cpp
    Source/Runtime/Renderer/Public/CpuImagePresenter.h Source/Runtime/Renderer/Private/CpuImagePresenter.cpp)
target_link_libraries(MooRenderer PRIVATE MooRHI glad::glad)

moo_module(MooApplicationCore Source/Runtime/ApplicationCore
    Source/Runtime/ApplicationCore/Public/Platform.h Source/Runtime/ApplicationCore/Private/Platform.cpp
    Source/Runtime/ApplicationCore/Public/AppWindow.h Source/Runtime/ApplicationCore/Private/AppWindow.cpp
    Source/Runtime/ApplicationCore/Public/Input.h Source/Runtime/ApplicationCore/Private/Input.cpp)
target_link_libraries(MooApplicationCore PRIVATE glfw glad::glad imgui::imgui)
if(WIN32)
    target_link_libraries(MooApplicationCore PRIVATE Xinput winmm)
endif()

# 纯数据动画算法：不依赖存量头或GL。
moo_module(MooAnimation Source/Runtime/Animation
    Source/Runtime/Animation/Public/FootContacts.h Source/Runtime/Animation/Private/FootContacts.cpp
    Source/Runtime/Animation/Public/TwoBoneIk.h Source/Runtime/Animation/Private/TwoBoneIk.cpp)

moo_module(MooSandbox Source/Sandbox
    Source/Sandbox/Public/SceneTypes.h Source/Sandbox/Public/SandboxScene.h
    Source/Sandbox/Private/SandboxScene.cpp Source/Sandbox/Private/SceneRender.inl Source/Sandbox/Private/MotionMatching.inl
    Source/Sandbox/Private/FootIk.inl)
set_source_files_properties(Source/Sandbox/Private/SceneRender.inl Source/Sandbox/Private/MotionMatching.inl
    Source/Sandbox/Private/FootIk.inl PROPERTIES HEADER_FILE_ONLY TRUE)
target_include_directories(MooSandbox PRIVATE Source/Legacy)
target_link_libraries(MooSandbox PUBLIC MooRenderer PRIVATE MooAnimation ${OpenCV_LIBS} OpenMP::OpenMP_CXX)
target_compile_definitions(MooSandbox PRIVATE MOOLAB_RESOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/resources")
moo_module(MooEditor Source/Editor Source/Editor/Public/EditorUi.h Source/Editor/Private/EditorUi.cpp)
target_include_directories(MooEditor PUBLIC Source/Sandbox/Public)
target_link_libraries(MooEditor PRIVATE imgui::imgui)
add_executable(MooLab Source/Launch/Private/Main.cpp Source/Launch/Private/Application.h Source/Launch/Private/Application.cpp)
target_link_libraries(MooLab PRIVATE MooApplicationCore MooSandbox MooRenderer MooEditor ${OpenCV_LIBS})

set(MOOLAB_SHADERS
    Shaders/mesh.vert Shaders/light.frag Shaders/pbr.frag Shaders/point.vert
    Shaders/point.frag Shaders/fullscreen.vert Shaders/present.frag Shaders/pbr_common.glsl Shaders/ground.frag
    Shaders/pbr_probe.frag)
target_sources(MooLab PRIVATE ${MOOLAB_SHADERS})
source_group("Shaders" FILES ${MOOLAB_SHADERS})
set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT MooLab)
set_target_properties(MooLab PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/resources")

add_executable(MooImageDiff Tools/ImageDiff/ImageDiff.cpp)
target_link_libraries(MooImageDiff PRIVATE MooCore ${OpenCV_LIBS})
add_executable(SceneProbe EXCLUDE_FROM_ALL Tools/Verify/SceneProbe.cpp)
target_link_libraries(SceneProbe PRIVATE MooSandbox MooAnimation MooApplicationCore ${OpenCV_LIBS})
add_executable(FootContactBake EXCLUDE_FROM_ALL Tools/FootContacts/FootContactBake.cpp)
target_link_libraries(FootContactBake PRIVATE MooSandbox)
add_executable(gl_probe EXCLUDE_FROM_ALL Tools/Verify/gl_probe.cpp)
target_link_libraries(gl_probe PRIVATE MooRenderer MooRHI MooApplicationCore glad::glad)
add_executable(MaterialProbe EXCLUDE_FROM_ALL Tools/Verify/MaterialProbe.cpp)
target_link_libraries(MaterialProbe PRIVATE MooRenderer MooRHI MooApplicationCore glad::glad)

# 公共头文件名必须唯一；目录与构建清单的增加均为显式注册。
set(MOO_PUBLIC_HEADERS
    Source/Runtime/Core/Public/Core.h Source/Runtime/Core/Public/GlContext.h Source/Runtime/Core/Public/ViewportTexture.h
    Source/Runtime/RHI/Public/GlRhi.h Source/Runtime/Renderer/Public/RenderTypes.h Source/Runtime/Renderer/Public/PbrShading.h
    Source/Runtime/Renderer/Public/GpuRenderer.h Source/Runtime/Renderer/Public/CpuImagePresenter.h
    Source/Runtime/ApplicationCore/Public/Platform.h Source/Runtime/ApplicationCore/Public/AppWindow.h
    Source/Sandbox/Public/SceneTypes.h Source/Sandbox/Public/SandboxScene.h
    Source/Runtime/ApplicationCore/Public/Input.h Source/Editor/Public/EditorUi.h
    Source/Runtime/Animation/Public/FootContacts.h Source/Runtime/Animation/Public/TwoBoneIk.h)
set(MOO_PUBLIC_NAMES)
foreach(Header IN LISTS MOO_PUBLIC_HEADERS)
    get_filename_component(Name "${Header}" NAME)
    if(Name IN_LIST MOO_PUBLIC_NAMES)
        message(FATAL_ERROR "Duplicate public header basename: ${Name}")
    endif()
    list(APPEND MOO_PUBLIC_NAMES "${Name}")
endforeach()
