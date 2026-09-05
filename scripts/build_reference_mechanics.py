#!/usr/bin/env python3
"""Build a pinned, repository-local FEBio reference; no system installation."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
COMMIT = "32ae206ff4881dfb54f62296cd1558e58ed9fcc6"
SOURCE = ROOT / "data/raw/mechanics/FEBio-4.13"
OUT = ROOT / "data/derived/mechanics-reference"
BUILD = OUT / "build"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    records = []
    report = {"backend": "FEBio", "version": "4.13", "commit": COMMIT,
              "source_url": "https://github.com/febiosoftware/FEBio.git",
              "machine": platform.machine(), "platform": platform.platform(),
              "commands": records, "status": "building"}
    def run(argv, name):
        started = time.time()
        with (OUT / name).open("w") as log:
            result = subprocess.run(argv, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
        records.append({"argv": argv, "returncode": result.returncode,
                        "elapsed_seconds": time.time() - started, "log": name})
        if result.returncode:
            raise RuntimeError(f"Command failed; see {OUT / name}")
    try:
        if not SOURCE.exists():
            SOURCE.parent.mkdir(parents=True, exist_ok=True)
            run(["git", "clone", "--depth", "1", "--branch", "v4.13", report["source_url"], str(SOURCE)], "clone.log")
        actual = subprocess.check_output(["git", "-C", str(SOURCE), "rev-parse", "HEAD"], text=True).strip()
        if actual != COMMIT:
            raise RuntimeError(f"Source revision mismatch: {actual}")
        if subprocess.check_output(["git", "-C", str(SOURCE), "status", "--porcelain"], text=True).strip():
            raise RuntimeError("Reference source must be unmodified")
        report["source_build_sha256"] = sha256(SOURCE / "CMakeLists.txt")
        archive = SOURCE.parent / 'FEBio-4.13.tar'
        run(['git', '-C', str(SOURCE), 'archive', '--format=tar', '-o', str(archive), COMMIT], 'source-archive.log')
        report['source_archive_sha256'] = sha256(archive)
        flags = [f"-DUSE_{name}=OFF" for name in ("MKL", "HYPRE", "SUPERLU_MT", "MMG", "LEVMAR", "NLOPT", "PDL", "FFTW", "ZLIB", "STATIC_STDLIBS")]
        # Upstream enables OpenMP compilation but its no-MKL Linux path omits
        # the runtime linkage. Supply driver flags without modifying sources.
        run(["cmake", "-S", str(SOURCE), "-B", str(BUILD), "-DCMAKE_BUILD_TYPE=Release",
             "-DCMAKE_EXE_LINKER_FLAGS=-fopenmp", "-DCMAKE_SHARED_LINKER_FLAGS=-fopenmp", *flags], "configure.log")
        run(["cmake", "--build", str(BUILD), "--parallel", str(args.jobs)], "build.log")
        binary = BUILD / "bin/febio4"
        report.update(status="built_not_benchmarked", executable=str(binary.relative_to(ROOT)), executable_sha256=sha256(binary))
        run(["ldd", str(binary)], "linked-libraries.log")
        report["shared_library_sha256"] = {str(p.relative_to(ROOT)): sha256(p) for p in sorted((BUILD / "lib").glob("*.so"))}
    except Exception as error:
        report.update(status="failed", error=str(error))
        raise
    finally:
        (OUT / "build-manifest.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
