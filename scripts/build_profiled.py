#!/usr/bin/env python3
"""Build the retained Rust code with the saved, workload-specific LLVM profile."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile', type=Path,
                        default=ROOT/'profiles/mumento2-10000-llvm22.profdata')
    parser.add_argument('--metadata', type=Path,
                        default=ROOT/'profiles/mumento2-10000-llvm22.json')
    parser.add_argument('-j', '--jobs', type=int, default=4)
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error('jobs must be positive')
    metadata = json.loads(args.metadata.read_text())
    profile = args.profile.resolve(strict=True)
    digest = hashlib.sha256(profile.read_bytes()).hexdigest()
    if digest != metadata['profile_sha256']:
        parser.error('profile digest does not match metadata')
    for name, expected in metadata['source_sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != expected:
            parser.error(f'{name} changed; retrain a profile for this source version')
    compiler = subprocess.check_output(['rustc', '-vV'], text=True)
    if compiler != metadata['rustc']:
        parser.error('saved profile requires the recorded compiler; retrain for this compiler')
    host = next(line.split(': ', 1)[1] for line in compiler.splitlines() if line.startswith('host: '))
    target = ROOT/'target/profiled'
    profiles = target/'profiles'
    profiles.mkdir(parents=True, exist_ok=True)
    frozen = profiles/(digest+'.profdata')
    shutil.copyfile(profile, frozen)
    env = os.environ.copy()
    env.pop('CARGO_ENCODED_RUSTFLAGS', None)
    env['RUSTFLAGS'] = f'-Ctarget-cpu=native -Cprofile-use={frozen} -Cllvm-args=-pgo-warn-missing-function'
    subprocess.run(['cargo', 'build', '--locked', '--release', '--target', host,
                    '--target-dir', str(target), '-j', str(args.jobs)],
                   cwd=ROOT, env=env, check=True)
    print(target/host/'release/mumento2')


if __name__ == '__main__':
    main()
