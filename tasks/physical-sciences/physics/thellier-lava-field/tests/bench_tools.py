"""Small root helpers for the bench: seal, restage, cpus, file, scrub (one per sub-command)."""
import hashlib
import json
import os
import shutil
import stat
import sys
from pathlib import Path

CAP = 65536
SOLUTION = Path("/app/solution")


def digests():
    return sorted(f"{hashlib.sha256(p.read_bytes()).hexdigest()} {p}" for p in Path("/tests").rglob("*")
                  if p.is_file() and "__pycache__" not in p.parts)


def seal(out):
    Path(out).write_text("\n".join(digests()))
    os.chmod(out, 0o600)


def restage(shelf):
    """rebuilds /app/solution in place from its own regular files and folders, root-owned and read-only."""
    shelf = Path(shelf)
    shelf.mkdir(mode=0o700)
    keep = []
    for base, folders, files in os.walk(SOLUTION):
        rel = Path(base).relative_to(SOLUTION)
        folders[:] = [f for f in folders if not (Path(base) / f).is_symlink()]
        keep += [(rel / f, True) for f in folders]
        keep += [(rel / f, False) for f in files if stat.S_ISREG(os.lstat(Path(base) / f).st_mode)]
    for rel, folder in keep:
        (shelf / rel).mkdir() if folder else shutil.copyfile(SOLUTION / rel, shelf / rel)
    for item in list(SOLUTION.iterdir()):
        shutil.rmtree(item) if item.is_dir() and not item.is_symlink() else item.unlink()
    for rel, folder in keep:
        (SOLUTION / rel).mkdir() if folder else shutil.copyfile(shelf / rel, SOLUTION / rel)
    for rel, folder in [(Path("."), True)] + keep:
        os.chown(SOLUTION / rel, 0, 0)
        os.chmod(SOLUTION / rel, 0o555 if folder else 0o444)
    shutil.rmtree(shelf)


def cpus():
    print(",".join(str(c) for c in sorted(os.sched_getaffinity(0))[:2]))


def file_answer(path, target, lava, code, left):
    outcome = "absent"
    try:
        fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    except OSError:
        fd = None
    if fd is not None:
        with os.fdopen(fd, "rb") as src:
            info = os.fstat(src.fileno())
            body = src.read(CAP + 1) if stat.S_ISREG(info.st_mode) else None
        if body is None or len(body) > CAP:
            outcome = "refused"
        else:
            Path(target).write_bytes(body)
            outcome = "filed"
    print(json.dumps({"lava": lava, "exit": code, "swept": left == "0", "answer": outcome}))


def scrub():
    for top in ("/tmp", "/var/tmp", "/dev/shm"):
        for item in Path(top).glob("*") if os.path.isdir(top) else ():
            info = item.lstat()
            if info.st_uid == 65534:
                shutil.rmtree(item, ignore_errors=True) if stat.S_ISDIR(info.st_mode) else item.unlink(missing_ok=True)


if __name__ == "__main__":
    {"seal": seal, "restage": restage, "cpus": cpus, "file": file_answer, "scrub": scrub}[sys.argv[1]](*sys.argv[2:])
