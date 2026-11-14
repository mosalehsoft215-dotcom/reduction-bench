"""drop_exec.py <cpu list> <home> <program> <casting copy> <target>: started by lift_spawn.py as root inside a new
session with its streams already on /dev/null; clears the supplementary groups, becomes gid/uid 65534, sets
no_new_privs, pins itself to the given CPUs, enters <home> and replaces itself with python3 running <program>, keeping
the environment it was given (PATH and HOME only).
"""
import ctypes, os, sys

NOBODY = 65534


def drop(cpus, home):
    os.setgroups([])
    os.setgid(NOBODY)
    os.setuid(NOBODY)
    if ctypes.CDLL(None).prctl(38, 1, 0, 0, 0) != 0:
        raise OSError("no_new_privs refused")
    os.sched_setaffinity(0, cpus)
    os.chdir(home)


if __name__ == "__main__":
    cpu_text, home, program, folder, target = sys.argv[1:6]
    try:
        drop({int(c) for c in cpu_text.split(",")}, home)
        os.execve("/usr/local/bin/python3", ["python3", program, folder, target], dict(os.environ))
    finally:
        os._exit(126)
