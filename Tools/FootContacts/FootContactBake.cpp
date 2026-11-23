// 由resources/long_motion.bin重算逐脚触地标签并写出foot_contacts.bin。
// 用法：FootContactBake [输出路径]，缺省写入资源目录。只新增标签文件，不改其他资源。
#include "SandboxScene.h"
#include <iostream>
#include <string>

int main(int Argc, char** Argv)
{
    const std::string Output = Argc > 1 ? Argv[1] : GetSandboxResourcePath("foot_contacts.bin");
    std::string Report;
    if(!BakeFootContacts(Output, Report))
    {
        std::cerr << "FootContactBake failed: " << Report << '\n';
        return 1;
    }
    std::cout << "wrote " << Output << '\n' << Report << '\n';
    return 0;
}
