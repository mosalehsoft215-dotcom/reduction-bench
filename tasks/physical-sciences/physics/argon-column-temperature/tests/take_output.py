"""Copy a program's output into the ledger only if it is a regular file of at most 65536 bytes, opened without
following links and without blocking (so a FIFO or a symlink is never read)."""
import os
import stat
import sys

src, dst = sys.argv[1], sys.argv[2]
try:
    fd = os.open(src, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
except OSError:
    sys.exit(0)
with os.fdopen(fd, "rb") as handle:
    info = os.fstat(handle.fileno())
    if stat.S_ISREG(info.st_mode) and info.st_size <= 65536:
        body = handle.read(65537)
        if len(body) <= 65536:
            with open(dst, "wb") as out:
                out.write(body)
