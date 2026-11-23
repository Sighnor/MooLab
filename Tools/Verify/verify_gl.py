import ctypes
from ctypes import wintypes
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

from PIL import Image

REPOSITORY = Path(__file__).resolve().parents[2]
BUILD_ARG = next((value.split('=', 1)[1] for value in sys.argv if value.startswith('--build-dir=')), 'build/vs2026')
ROOT = (REPOSITORY / BUILD_ARG).resolve()
EVIDENCE_ARG = next((value.split('=', 1)[1] for value in sys.argv if value.startswith('--evidence-dir=')), 'foot_ik_v1')
# 历史地面证据只读；新证据只允许位于明确的foot_ik版本目录，拒绝绝对路径和目录穿越。
if not re.fullmatch(r'foot_ik_v[1-9]\d*', EVIDENCE_ARG):
    raise ValueError('--evidence-dir只能是foot_ik_vN；ground_pbr_v1及历史证据只读')
EVIDENCE = ROOT / EVIDENCE_ARG
GROUND_BASELINE = ROOT / 'ground_pbr_v1/baseline'
CANDIDATE = EVIDENCE / 'candidate'
BASELINE = EVIDENCE / 'baseline'
EXE = ROOT / 'Release/MooLab.exe'
DIFF = ROOT / 'Release/MooImageDiff.exe'


def run_case(name, args, expected=0, exe=EXE, env=None, timeout=180):
    with (CANDIDATE / f'{name}.log').open('wb') as out, (CANDIDATE / f'{name}.err').open('wb') as err:
        process = subprocess.run([str(exe), *args], stdout=out, stderr=err, env=env, timeout=timeout)
        assert process.returncode == expected, (name, process.returncode)
    print(f'{name}: exit={process.returncode}', flush=True)


def capture():
    check_legacy()
    inputs = input_hashes()
    env = dict(os.environ, OMP_WAIT_POLICY='PASSIVE')
    for prefix, mode, config, hidden, threads in [
        ('cpu', 'cpu', 'Release', True, None),
        ('gpu', 'gpu', 'Release', True, None),
        ('repeat_gpu', 'gpu', 'Release', True, None),
        ('visible_gpu', 'gpu', 'Release', False, None),
        ('debug_gpu', 'gpu', 'Debug', True, None),
        ('debug_cpu', 'cpu', 'Debug', True, None),
        ('single_cpu', 'cpu', 'Release', True, '1'),
    ]:
        args = ['--renderer', mode, '--frames', '301', '--unlimited-fps']
        if hidden:
            args.append('--hidden')
        for frame in [60, 300]:
            args += ['--capture-frame', str(frame), str(CANDIDATE / f'{prefix}_{frame:03}.png')]
        case_env = dict(env)
        if threads:
            case_env['OMP_NUM_THREADS'] = threads
        run_case(f'capture_{prefix}', args, exe=ROOT / config / 'MooLab.exe', env=case_env)
        assert (CANDIDATE / f'capture_{prefix}.err').stat().st_size == 0, prefix
    files = {f'{prefix}_{frame:03}.png': sha(CANDIDATE / f'{prefix}_{frame:03}.png')
             for prefix in ['cpu', 'gpu', 'repeat_gpu', 'visible_gpu', 'debug_gpu', 'debug_cpu', 'single_cpu'] for frame in [60, 300]}
    assert inputs == input_hashes(), 'inputs changed during capture'
    (CANDIDATE / 'capture.json').write_text(json.dumps(dict(inputs=inputs, images=files), indent=2), encoding='utf-8')
    verify_candidate()


def probe_config(config):
    exe = ROOT / config / 'gl_probe.exe'
    run_case(f'gl_probe_{config}', [], exe=exe)
    assert (CANDIDATE / f'gl_probe_{config}.err').stat().st_size == 0
    for name, argument, diagnostic in [
        ('probe_bad_context', 'bad-context', b'context'),
        ('probe_bad_update', 'bad-update', b'vertex count'),
        ('probe_bad_index', 'bad-index', b'out of range'),
        ('probe_bad_present', 'bad-present', b'does not match presenter'),
        ('ground_outside', 'bad-ground-outside', b'outside BeginFrame'),
        ('ground_twice', 'bad-ground-twice', b'at most once'),
        ('ground_after_mesh', 'bad-ground-after-mesh', b'first draw'),
        ('ground_after_points', 'bad-ground-after-points', b'first draw'),
        ('ground_size', 'bad-ground-size', b'do not match target'),
    ]:
        name += '_' + config
        with (CANDIDATE / f'{name}.log').open('wb') as out, (CANDIDATE / f'{name}.err').open('wb') as err:
            result = subprocess.run([str(exe), argument], stdout=out, stderr=err, timeout=20)
        assert result.returncode != 0, name
        assert diagnostic in (CANDIDATE / f'{name}.err').read_bytes(), name
        print(f'{name}: expected failure exit={result.returncode}', flush=True)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def read_pre_implementation():
    data = (EVIDENCE / 'pre_implementation_sha256.txt').read_bytes()
    if data[:2] in (b'\xff\xfe', b'\xfe\xff'):
        return data.decode('utf-16')
    try:
        return data.decode('utf-8-sig')
    except UnicodeDecodeError:
        # 实施前清单由Windows PowerShell生成（CP936）；只读兼容，不改写旧证据
        return data.decode('cp936')


def check_legacy():
    entries = {}
    for line in read_pre_implementation().splitlines():
        if re.match(r'^[0-9A-Fa-f]{64} ', line):
            digest, name = line.split(' ', 1)
            entries[name] = digest.upper()
    assert entries and any(name.startswith('editor_baseline/') for name in entries), 'missing old evidence manifest'
    for name, digest in entries.items():
        assert sha(ROOT / name) == digest, ('old evidence changed', name)
    old_files = {p.relative_to(ROOT).as_posix() for p in (ROOT / 'editor_baseline').rglob('*') if p.is_file()}
    assert old_files == {name for name in entries if name.startswith('editor_baseline/')}, 'old baseline file set changed'
    ground_files = {p.relative_to(ROOT).as_posix() for p in GROUND_BASELINE.rglob('*') if p.is_file()}
    assert ground_files and ground_files == {name for name in entries if name.startswith('ground_pbr_v1/baseline/')}, 'ground baseline file set changed'
    print(f'old evidence SHA256: PASS ({len(entries)} files, ground baseline read-only)', flush=True)


def input_hashes():
    # 候选必须由当前源码/着色器与所选构建的可执行文件生成，防止旧候选被意外冻结
    files = [REPOSITORY / p for p in ['CMakeLists.txt', 'CMakePresets.json', 'vcpkg.json']]
    for folder in ['Source', 'Shaders', 'Tools']:
        files += [p for p in (REPOSITORY / folder).rglob('*') if p.is_file() and p.suffix in ['.h', '.hpp', '.inl', '.cpp', '.py', '.frag', '.vert', '.glsl']]
    files += [REPOSITORY / 'resources' / name for name in ['character.bin', 'database.bin', 'long_motion.bin', 'skybox.png', 'foot_contacts.bin']]
    hashes = {p.relative_to(REPOSITORY).as_posix(): sha(p) for p in files}
    for config in ['Debug', 'Release']:
        for name in ['MooLab', 'MooImageDiff', 'SceneProbe', 'gl_probe', 'MaterialProbe', 'FootContactBake']:
            path = ROOT / config / (name + '.exe')
            hashes[f'build/{config}/{name}.exe'] = sha(path)
    return hashes


def diff_images(a, b, expected=0):
    assert subprocess.run([str(DIFF), str(CANDIDATE / a), str(CANDIDATE / b)], timeout=30).returncode == expected, (a, b)


def check_corners(name, ground):
    with Image.open(CANDIDATE / name) as image:
        image = image.convert('RGB')
        assert image.size == (720, 720)
        assert all(image.getpixel(p) == (100, 100, 200) for p in [(0, 0), (719, 0)]), (name, 'upper background')
        assert all((image.getpixel(p) != (100, 100, 200)) == ground for p in [(0, 719), (719, 719)]), (name, 'lower ground partition')


def check_green(a, b):
    with Image.open(CANDIDATE / a) as gpu, Image.open(CANDIDATE / b) as cpu:
        gpu, cpu = gpu.convert('RGB'), cpu.convert('RGB')
        assert gpu.size == cpu.size == (720, 720)
        green_gpu = [p == (0, 255, 0) for p in gpu.get_flattened_data()]
        green_cpu = [p == (0, 255, 0) for p in cpu.get_flattened_data()]
        assert green_gpu == green_cpu, (a, b, 'contact mask')
        print(f'{a}: green-count={sum(green_gpu)} mismatch=0', flush=True)


def check_images():
    receipt = json.loads((CANDIDATE / 'capture.json').read_text(encoding='utf-8'))
    assert receipt['inputs'] == input_hashes(), 'source/binaries changed: run capture again'
    required = {f'{prefix}_{frame}.png' for prefix in ['cpu', 'gpu', 'repeat_gpu', 'visible_gpu', 'debug_gpu', 'debug_cpu', 'single_cpu'] for frame in ['060', '300']}
    assert set(receipt['images']) == required
    for name, digest in receipt['images'].items():
        assert sha(CANDIDATE / name) == digest, ('capture hash', name)
        assert digest == sha(GROUND_BASELINE / name), ('Physics image changed from ground baseline', name)
    for frame in ['060', '300']:
        for base, variants in [('cpu', ['debug_cpu', 'single_cpu']), ('gpu', ['repeat_gpu', 'visible_gpu', 'debug_gpu'])]:
            for variant in variants:
                assert sha(CANDIDATE / f'{variant}_{frame}.png') == sha(CANDIDATE / f'{base}_{frame}.png'), variant
        diff_images(f'gpu_{frame}.png', f'cpu_{frame}.png')
        check_green(f'gpu_{frame}.png', f'cpu_{frame}.png')
        for mode in ['gpu', 'cpu']:
            check_corners(f'{mode}_{frame}.png', True)
        print(f'frame {frame}: 7-config determinism / ground corners: PASS', flush=True)
    diff_images('gpu_060.png', 'cpu_300.png', 1)
    print('cross-frame negative control: PASS', flush=True)


def cli_negatives():
    for name, args in [('invalid', ['--renderer', 'invalid']), ('incomplete_frames', ['--hidden', '--frames']),
                       ('incomplete_capture', ['--hidden', '--capture-frame', '5']),
                       ('invalid_scene', ['--hidden', '--scene', 'invalid']), ('incomplete_scene', ['--hidden', '--scene'])]:
        with (CANDIDATE / f'{name}.log').open('wb') as out, (CANDIDATE / f'{name}.err').open('wb') as err:
            result = subprocess.run([str(EXE), *args], stdout=out, stderr=err, timeout=20)
        assert result.returncode != 0, name
        print(f'{name} exit={result.returncode}', flush=True)
    run_case('unreached', ['--hidden', '--frames', '0', '--capture-frame', '0', str(CANDIDATE / 'never.png')], 1)
    run_case('write_failed', ['--hidden', '--frames', '1', '--capture-frame', '0', str(CANDIDATE / 'nonexistent-directory/capture.png')], 1)
    run_case('write_exception', ['--hidden', '--frames', '1', '--capture-frame', '0', str(CANDIDATE / 'capture.unsupported-format')], 1)
    print('CLI negative cases: PASS', flush=True)


def probe():
    for config in ['Release', 'Debug']:
        probe_config(config)


def brdf(values, mutant=None):
    # Python double独立求值，不调用生产数值核或从CPU/GPU结果反推
    albedo, r, m = values[:3], values[3], values[4]
    nv, nl, nh, vh = values[5:]
    if nv <= 0 or nl <= 0:
        return [0.] * 6
    r, m = min(1., max(.05, r)), min(1., max(0., m))
    nv, nl, nh, vh = [min(1., max(0., x)) for x in [nv, nl, nh, vh]]
    alpha2 = r ** 4
    d = alpha2 / (math.pi * ((1 - nh * nh) + nh * nh * alpha2) ** 2)
    k = (r + 1) ** 2 / 8
    g = nv / (nv * (1 - k) + k) * nl / (nl * (1 - k) + k)
    f0 = [(.04 * (1 - m) + m * a) if mutant != 'f0' else .04 for a in albedo]
    cosine = nv if mutant == 'fresnel' else vh
    fresnel = [f + (1 - f) * (1 - cosine) ** 5 for f in f0]
    diffuse = [(1 if mutant == 'metal-diffuse' else 1 - m) * (1 - f) * a / math.pi for a, f in zip(albedo, fresnel)]
    specular = [d * g * f / max(4 * nv * nl, .001) for f in fresnel]
    return diffuse + specular


def close(actual, reference, absolute, relative):
    return all(math.isfinite(a) and abs(a - b) <= absolute + relative * abs(b) for a, b in zip(actual, reference))


def material():
    for config in ['Release', 'Debug']:
        name = 'material_' + config
        run_case(name, [], exe=ROOT / config / 'MaterialProbe.exe')
        assert (CANDIDATE / f'{name}.err').stat().st_size == 0
        samples = {}
        text = (CANDIDATE / f'{name}.log').read_text(encoding='utf-8')
        assert 'material-probe samples=46: PASS' in text
        for line in text.splitlines():
            if not line.startswith('SAMPLE '):
                continue
            fields = line.split()
            assert len(fields) == 23 and fields[1] not in samples, line
            numbers = [float(x) for x in fields[2:]]
            values, cpu, gpu = numbers[:9], numbers[9:15], numbers[15:21]
            reference = brdf(values)
            assert close(cpu, reference, 1e-5, 1e-5), (fields[1], 'CPU vs Python', cpu, reference)
            assert close(gpu, reference, 2e-5, 5e-4), (fields[1], 'GPU vs Python', gpu, reference)
            if values[4] == 1:
                assert cpu[:3] == gpu[:3] == [0., 0., 0.], (fields[1], 'metal diffuse not exactly zero')
            samples[fields[1]] = values
        assert len(samples) == 46
        for m in [0, 5, 10]:
            for r in [5, 50, 100]:
                for direction in range(4):
                    assert f'grid_m{m}_r{r}_d{direction}' in samples
        for tag, anchor in [('anchor_dielectric', .3565071), ('anchor_metal', 1.2732395), ('anchor_half', .4018662), ('anchor_fresnel', .0042485)]:
            reference = brdf(samples[tag])
            assert abs(reference[0] + reference[3] - anchor) < 1e-6, tag
        for mutant, tag in [('f0', 'anchor_metal'), ('fresnel', 'anchor_fresnel'), ('metal-diffuse', 'bound_r1_m1')]:
            values = samples[tag]
            assert not close(brdf(values, mutant), brdf(values), 2e-5, 5e-4), (mutant, 'negative control not rejected')
        print(f'{config} material: 46 samples / independent double / anchors / 3 mutant controls: PASS', flush=True)


def scene(use_cache=False):
    inputs = input_hashes()
    names = ['scene_release', 'scene_debug', 'scene_single']
    suffixes = ['gpu.png', 'cpu.png', 'hidden_gpu.png', 'hidden_cpu.png', 'mm_gpu_060.png', 'mm_cpu_060.png', 'mm_gpu_300.png', 'mm_cpu_300.png']
    required = {f'{name}_{suffix}' for name in names for suffix in suffixes} | {f'{name}.{ext}' for name in names for ext in ['log', 'err']}
    receipt_path = CANDIDATE / 'scene_checks.json'
    if use_cache and receipt_path.exists():
        receipt = json.loads(receipt_path.read_text(encoding='utf-8'))
        if receipt['inputs'] == inputs and set(receipt['files']) == required and all((CANDIDATE / name).is_file() and sha(CANDIDATE / name) == digest for name, digest in receipt['files'].items()):
            print('scene: prior full run source/binary/evidence SHA256 revalidated: PASS', flush=True)
            return
    evidence_files = []
    reference_digest = None
    for name, config, threads in [('scene_release', 'Release', None), ('scene_debug', 'Debug', None), ('scene_single', 'Release', '1')]:
        env = dict(os.environ, OMP_WAIT_POLICY='PASSIVE')
        if threads:
            env['OMP_NUM_THREADS'] = threads
        # Debug执行6万步随机复跑不优化，给予独立的工具看门狗；不改变任何运动预算
        run_case(name, [str(CANDIDATE / name)], exe=ROOT / config / 'SceneProbe.exe', env=env, timeout=900)
        # 缺失/损坏标签负对照会各输出一条预期MM不可用诊断；其余stderr仍拒绝。
        errors = (CANDIDATE / f'{name}.err').read_text(encoding='utf-8').splitlines()
        assert len(errors) == 2 and all(line.startswith('Motion Matching unavailable: ') and '逐脚触地标签无效' in line for line in errors), (name, errors)
        text = (CANDIDATE / f'{name}.log').read_text(encoding='utf-8')
        assert 'hard checks=PASS budgetFailures=0' in text
        old = re.search(r'motion digest physics600=([0-9a-f]+) matching3000=([0-9a-f]+)', read_pre_implementation())
        current = re.search(r'motion digest default physics600=([0-9a-f]+) matching3000=([0-9a-f]+) ground-edited physics600=([0-9a-f]+) matching3000=([0-9a-f]+)', text)
        assert old and current and current[1] == current[3] == old[1] and current[2] == current[4], 'Physics or ground-edit motion digest changed'
        digest = dict(physics600=current[1], matching3000=current[2])
        if reference_digest is None:
            reference_digest = digest
        assert digest == reference_digest, (name, 'Debug/Release/single-thread motion digest mismatch')
        for suffix in ['gpu.png', 'cpu.png', 'hidden_gpu.png', 'hidden_cpu.png']:
            assert sha(CANDIDATE / f'{name}_{suffix}') == sha(GROUND_BASELINE / f'{name}_{suffix}'), (name, suffix, 'Physics baseline mismatch')
        for mode in ['cpu', 'gpu']:
            check_corners(f'{name}_hidden_{mode}.png', False)
        for frame in ['060', '300']:
            diff_images(f'{name}_mm_gpu_{frame}.png', f'{name}_mm_cpu_{frame}.png')
            check_green(f'{name}_mm_gpu_{frame}.png', f'{name}_mm_cpu_{frame}.png')
        evidence_files += [CANDIDATE / f'{name}.log', CANDIDATE / f'{name}.err']
        for suffix in suffixes:
            evidence_files.append(CANDIDATE / f'{name}_{suffix}')
            assert sha(CANDIDATE / f'{name}_{suffix}') == sha(CANDIDATE / f'scene_release_{suffix}'), (name, suffix, 'scene determinism')
    assert inputs == input_hashes(), 'inputs changed during scene checks'
    receipt_path.write_text(json.dumps(dict(inputs=inputs, digests=reference_digest, files={p.name: sha(p) for p in evidence_files}), indent=2), encoding='utf-8')
    print('scene geometry / settings / motion / budgets / MM images / no-ground corners: PASS', flush=True)


def verify_candidate():
    check_legacy()
    check_images()
    material()
    probe()
    scene(use_cache=True)
    cli_negatives()
    check_legacy()
    print('candidate binary checks: PASS; no baseline needed', flush=True)


def check_baseline():
    assert BASELINE.is_dir(), 'new baseline missing: capture, verify-candidate, then explicit freeze-baseline'
    manifest = json.loads((BASELINE / 'sha256.json').read_text(encoding='utf-8'))
    actual = {p.relative_to(BASELINE).as_posix() for p in BASELINE.rglob('*') if p.is_file() and p.name != 'sha256.json'}
    assert set(manifest) == actual and 'cpu_060.png' in manifest and 'gpu_300.png' in manifest, 'baseline file set wrong'
    for name, digest in manifest.items():
        assert sha(BASELINE / name) == digest, ('frozen baseline modified', name)
    return manifest


def freeze_baseline():
    assert not BASELINE.exists(), 'baseline already exists; overwriting is forbidden'
    verify_candidate()
    capture_receipt = json.loads((CANDIDATE / 'capture.json').read_text(encoding='utf-8'))
    scene_receipt = json.loads((CANDIDATE / 'scene_checks.json').read_text(encoding='utf-8'))
    files = set(capture_receipt['images']) | set(scene_receipt['files']) | {'capture.json', 'scene_checks.json'}
    cases = [f'capture_{p}' for p in ['cpu', 'gpu', 'repeat_gpu', 'visible_gpu', 'debug_gpu', 'debug_cpu', 'single_cpu']]
    for config in ['Release', 'Debug']:
        cases += [f'{name}_{config}' for name in ['material', 'gl_probe', 'probe_bad_context', 'probe_bad_update', 'probe_bad_index',
                                                'probe_bad_present', 'ground_outside', 'ground_twice', 'ground_after_mesh', 'ground_after_points', 'ground_size']]
    cases += ['invalid', 'incomplete_frames', 'incomplete_capture', 'invalid_scene', 'incomplete_scene', 'unreached', 'write_failed', 'write_exception']
    files |= {f'{name}.{ext}' for name in cases for ext in ['log', 'err']}
    perf_receipt_path = CANDIDATE / 'perf_checks.json'
    if perf_receipt_path.exists():
        perf_receipt = json.loads(perf_receipt_path.read_text(encoding='utf-8'))
        assert perf_receipt['inputs'] == input_hashes(), 'performance evidence stale'
        for name, digest in perf_receipt['files'].items():
            assert sha(CANDIDATE / name) == digest, ('performance evidence changed', name)
        files |= set(perf_receipt['files']) | {'perf_checks.json'}
    assert all((CANDIDATE / name).is_file() for name in files), 'missing verified evidence'
    # mkdir拒绝已存在目录，只复制通过验证的清单；早期调试/失败文件不成为基线证据
    BASELINE.mkdir()
    for name in sorted(files):
        shutil.copy2(CANDIDATE / name, BASELINE / name)
    manifest = {name: sha(BASELINE / name) for name in sorted(files)}
    (BASELINE / 'sha256.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    check_baseline()
    print(f'new baseline frozen: {BASELINE}', flush=True)


def verify():
    manifest = check_baseline()
    verify_candidate()
    for name, digest in manifest.items():
        if name.endswith('.png'):
            assert sha(CANDIDATE / name) == digest, ('new image differs from frozen baseline', name)
    print('new baseline full SHA256 / candidate comparison: PASS', flush=True)


def percentile(values, p):
    values = sorted(values)
    return values[math.ceil(len(values) * p) - 1]


def perf():
    inputs = input_hashes()
    kernel = ctypes.WinDLL('kernel32', use_last_error=True)
    kernel.GetProcessTimes.argtypes = [wintypes.HANDLE] + [ctypes.POINTER(wintypes.FILETIME)] * 4
    kernel.GetProcessTimes.restype = wintypes.BOOL
    env = dict(os.environ, OMP_WAIT_POLICY='PASSIVE')
    pattern = re.compile(r'Frame: (\d+) Cost: ([\d.]+)ms Interval: ([\d.]+)ms Gpu: ([\d.]+)ms')
    samples = []
    for mode in ['gpu', 'cpu']:
        for index in range(1, 4):
            name = f'perf_{mode}_{index}'
            with (CANDIDATE / f'{name}.log').open('wb') as out, (CANDIDATE / f'{name}.err').open('wb') as err:
                clock = time.perf_counter()
                process = subprocess.Popen([str(EXE), '--renderer', mode, '--hidden', '--frames', '600'], stdout=out, stderr=err, env=env)
                assert process.wait(timeout=180) == 0
                wall = time.perf_counter() - clock
                stamps = [wintypes.FILETIME() for _ in range(4)]
                assert kernel.GetProcessTimes(int(process._handle), *(ctypes.byref(stamp) for stamp in stamps))
                cpu_seconds = sum(((stamp.dwHighDateTime << 32) | stamp.dwLowDateTime) / 1e7 for stamp in stamps[2:])
            rows = [(int(m[1]), float(m[2]), float(m[3]), float(m[4])) for m in pattern.finditer((CANDIDATE / f'{name}.log').read_text(encoding='utf-8'))]
            assert len(rows) == 600
            stable = [row for row in rows if row[0] >= 10]
            sample = dict(mode=mode, run=index, frames=600, stable_frames=len(stable),
                          cost_p50=percentile([r[1] for r in stable], .5), cost_p95=percentile([r[1] for r in stable], .95),
                          cost_p99=percentile([r[1] for r in stable], .99), gpu_p50=percentile([r[3] for r in stable], .5),
                          gpu_p95=percentile([r[3] for r in stable], .95), interval_p50=percentile([r[2] for r in rows[1:]], .5),
                          interval_p99=percentile([r[2] for r in rows[1:]], .99), process_cpu_seconds=cpu_seconds,
                          process_wall_seconds=wall, process_equivalent_cores=cpu_seconds / wall)
            samples.append(sample)
            print(json.dumps(sample), flush=True)
    # 进程等效核包含启动加载；稳定阶段Cost/GPU排除前10帧，不可混淆两种口径。
    (CANDIDATE / 'perf_gl.json').write_text(json.dumps(samples, indent=2), encoding='utf-8')
    run_case('ground_perf', [str(CANDIDATE / 'ground_perf'), '--ground-perf'], exe=ROOT / 'Release/SceneProbe.exe', env=env, timeout=900)
    assert (CANDIDATE / 'ground_perf.err').stat().st_size == 0
    groups = {}
    for line in (CANDIDATE / 'ground_perf.log').read_text(encoding='utf-8').splitlines():
        if not line.startswith('GROUND_PERF '):
            continue
        _, mode, run, visible, frame, cost, gpu = line.split()
        groups.setdefault((mode, int(run), int(visible)), []).append((int(frame), float(cost), float(gpu)))
    assert len(groups) == 12
    ground_samples = []
    for (mode, run, visible), rows in sorted(groups.items()):
        assert [r[0] for r in rows] == list(range(600))
        stable = rows[10:]
        gpu_times = [r[2] for r in stable if r[2] > 0.]
        assert mode == 'cpu' or len(gpu_times) >= 580, 'GPU timestamp samples not ready'
        sample = dict(mode=mode, run=run, ground_visible=bool(visible), frames=600, stable_frames=590,
                      draw_p50=percentile([r[1] for r in stable], .5), draw_p95=percentile([r[1] for r in stable], .95),
                      gpu_p50=percentile(gpu_times, .5) if gpu_times else 0., gpu_p95=percentile(gpu_times, .95) if gpu_times else 0.,
                      gpu_ready_samples=len(gpu_times))
        ground_samples.append(sample)
        print(json.dumps(sample), flush=True)
    (CANDIDATE / 'ground_perf.json').write_text(json.dumps(ground_samples, indent=2), encoding='utf-8')
    assert inputs == input_hashes(), 'inputs changed during performance checks'
    names = [f'perf_{mode}_{run}.{ext}' for mode in ['gpu', 'cpu'] for run in [1, 2, 3] for ext in ['log', 'err']]
    names += ['perf_gl.json', 'ground_perf.json', 'ground_perf.log', 'ground_perf.err']
    (CANDIDATE / 'perf_checks.json').write_text(json.dumps(dict(inputs=inputs, files={name: sha(CANDIDATE / name) for name in names}), indent=2), encoding='utf-8')


if __name__ == '__main__':
    if sys.flags.optimize:
        raise RuntimeError('验证脚本禁止 -O：判据必须启用')
    CANDIDATE.mkdir(parents=True, exist_ok=True)
    command = next((arg for arg in sys.argv[1:] if not arg.startswith(('--build-dir=', '--evidence-dir='))), 'verify')
    {'history': check_legacy, 'perf': perf, 'capture': capture, 'probe': probe, 'material': material, 'scene': scene,
     'verify-candidate': verify_candidate, 'freeze-baseline': freeze_baseline, 'verify': verify}[command]()
