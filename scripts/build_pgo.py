#!/usr/bin/env python3
"""Build native, profile-guided Rust and GGCAT binaries using representative genomes."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]


def run(command, *, env=None):
    print('+ ' + shlex.join(map(str, command)), flush=True)
    subprocess.run(list(map(str, command)), cwd=ROOT, env=env, check=True)


def profile_environment(option):
    env = os.environ.copy()
    if 'CARGO_ENCODED_RUSTFLAGS' in env:
        flags = [flag for flag in env['CARGO_ENCODED_RUSTFLAGS'].split('\x1f') if flag]
    else:
        flags = shlex.split(env.get('RUSTFLAGS', '-C target-cpu=native'))
    flags.append('-C' + option)
    env['CARGO_ENCODED_RUSTFLAGS'] = '\x1f'.join(flags)
    # Use the separate profile directories embedded in each instrumented binary.
    env.pop('LLVM_PROFILE_FILE', None)
    return env


def build(target, directory, *, graph, option, jobs):
    profile = 'mumento' if graph else 'release'
    manifest = ROOT / ('ggcat/Cargo.toml' if graph else 'Cargo.toml')
    command = ['cargo', 'build', '--locked', '--manifest-path', manifest,
               '--profile', profile, '--target', target, '--target-dir', directory,
               '-j', str(jobs)]
    if graph:
        command += ['-p', 'ggcat_cmdline', '--bin', 'ggcat']
    else:
        command += ['--bin', 'mumento2']
    run(command, env=profile_environment(option))
    return directory / target / profile / ('ggcat' if graph else 'mumento2')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('genomes', type=Path, nargs='+', help='Representative training inputs')
    parser.add_argument('-j', '--threads', type=int, default=os.cpu_count() or 1)
    parser.add_argument('-k', type=int, default=31)
    parser.add_argument('-m', '--memory', type=float, default=1.0)
    parser.add_argument('--buckets-count-log', type=int, default=8)
    parser.add_argument('--ggcat-minimizer-length', type=int)
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'target/pgo')
    args = parser.parse_args()
    if args.threads < 1 or args.k < 3:
        parser.error('threads must be positive and k must be at least 3')
    inputs = [path.resolve(strict=True) for path in args.genomes]
    if not all(path.is_file() for path in inputs):
        parser.error('training inputs must be files')
    target = next(line.split(': ', 1)[1] for line in subprocess.check_output(
        ['rustc', '-vV'], text=True).splitlines() if line.startswith('host: '))
    sysroot = Path(subprocess.check_output(['rustc', '--print', 'sysroot'], text=True).strip())
    profdata = sysroot / 'lib/rustlib' / target / 'bin/llvm-profdata'
    if not profdata.is_file():
        parser.error('Install the matching profiling tools: rustup component add llvm-tools')
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    profiles = {}
    with tempfile.TemporaryDirectory(prefix='training-', dir=output) as temporary:
        training = Path(temporary)
        generated = {}
        for name, graph in [('ggcat', True), ('mumento2', False)]:
            raw = training / name
            raw.mkdir()
            generated[name] = build(target, output / 'generate' / name,
                                    graph=graph, option=f'profile-generate={raw}',
                                    jobs=args.threads)
            profiles[name] = raw
        command = [generated['mumento2'], '--ggcat', generated['ggcat'], '-k', str(args.k),
                   '-j', str(args.threads), '-m', str(args.memory),
                   '--buckets-count-log', str(args.buckets_count_log),
                   '-o', training / 'matches.tsv.zst', '--temp-dir', training / 'scratch']
        if args.ggcat_minimizer_length is not None:
            command += ['--ggcat-minimizer-length', str(args.ggcat_minimizer_length)]
        command += inputs
        training_env = os.environ.copy()
        training_env.pop('LLVM_PROFILE_FILE', None)
        run(command, env=training_env)
        for name, directory in profiles.items():
            raw = sorted(directory.glob('*.profraw'))
            if not raw:
                raise RuntimeError(f'No {name} execution profiles were generated')
            merged = output / f'{name}.profdata'
            run([profdata, 'merge', '-o', merged, *raw])
            # Cargo fingerprints flags rather than the contents of external
            # .profdata files. A content-specific path invalidates stale builds.
            digest = hashlib.sha256(merged.read_bytes()).hexdigest()[:16]
            profile = output / f'{name}-{digest}.profdata'
            merged.replace(profile)
            profiles[name] = profile
    binaries = {}
    for name, graph in [('ggcat', True), ('mumento2', False)]:
        compiled = build(target, output / 'use' / name, graph=graph,
                         option=f'profile-use={profiles[name]}', jobs=args.threads)
        binaries[name] = output / name
        shutil.copy2(compiled, binaries[name])
    report = dict(training_genomes=[str(path) for path in inputs], k=args.k,
                  threads=args.threads, target=target, rustc=subprocess.check_output(
                      ['rustc', '-vV'], text=True), binaries={k: str(v) for k, v in binaries.items()},
                  profiles={k: str(v) for k, v in profiles.items()},
                  elapsed_seconds=time.monotonic() - started)
    (output / 'build.json').write_text(json.dumps(report, indent=2) + '\n')
    print('\nPGO build complete. Run with:', flush=True)
    print(shlex.join([str(binaries['mumento2']), '--ggcat', str(binaries['ggcat']),
                      '-k', str(args.k), '-j', str(args.threads), '-o', 'matches.tsv.zst',
                      *map(str, inputs)]), flush=True)


if __name__ == '__main__':
    main()
