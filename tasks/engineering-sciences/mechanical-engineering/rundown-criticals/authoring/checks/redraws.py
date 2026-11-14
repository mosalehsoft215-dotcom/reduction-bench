"""Fresh instrument noise for all seven units, several times over, reduced by the reference and by the second
solver (the replica only on the bundle noise, level 0).

    python redraws.py <scratch folder> <csv out> <level> [level ...]
"""
import json, subprocess, sys
from pathlib import Path

TOP = Path(__file__).resolve().parents[2]
GENERATOR = TOP / "authoring" / "units" / "make_rundowns.py"
PROGRAMS = {"reference": TOP / "solution" / "criticals.py",
            "second_solver": TOP / "authoring" / "checks" / "second_solver.py",
            "replica_low": TOP / "authoring" / "checks" / "replica_low" / "criticals.py"}


def answer(program, unit_dir, slot):
    done = subprocess.run([sys.executable, str(program), str(unit_dir), str(slot)], capture_output=True)
    return json.loads(slot.read_text()) if done.returncode == 0 else {}


def main():
    scratch, out = Path(sys.argv[1]), Path(sys.argv[2])
    lines = ["level,program,unit,key,found,true,error_pct"]
    for level in ["0"] + sys.argv[3:]:
        folder = scratch / f"level_{level}"
        subprocess.run([sys.executable, str(GENERATOR), "--noise", level, "--out", str(folder)], check=True,
                       capture_output=True)
        truth = json.loads((folder / "truth.json").read_text())
        for program, path in PROGRAMS.items():
            if program == "replica_low" and level != "0":
                continue
            for unit in sorted(truth):
                got = answer(path, folder / unit, folder / "answer.json")
                for key, true in truth[unit].items():
                    found = got.get(key, float("nan"))
                    lines.append(f"{level},{program},{unit},{key},{found},{true},{100 * (found / true - 1):+.3f}")
        out.write_text("\n".join(lines) + "\n", newline="\n")
        print("level", level, "done", flush=True)


if __name__ == "__main__":
    main()
