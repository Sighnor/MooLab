@echo off
set "VCPKG_ROOT=C:\Program Files\Microsoft Visual Studio\18\Community\VC\vcpkg"
pushd "%~dp0..\.."
"%VCPKG_ROOT%\..\..\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe" --build --preset vs2026-release --target gl_probe
set "RESULT=%errorlevel%"
popd
exit /b %RESULT%
