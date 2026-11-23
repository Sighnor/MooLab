import subprocess
import sys
from pathlib import Path

# 启动最终Release进行人工/ComputerUse交互，并落盘退出码。
ROOT = Path(__file__).resolve().parents[2] / 'build/vs2026'
name = sys.argv[1]
with (ROOT / f'{name}.log').open('wb') as out, (ROOT / f'{name}.err').open('wb') as err:
    code = subprocess.run([str(ROOT / 'Release/MooLab.exe'), '--renderer', 'gpu', *sys.argv[2:]], stdout=out, stderr=err).returncode
(ROOT / f'{name}.exit').write_text(str(code), encoding='ascii')
