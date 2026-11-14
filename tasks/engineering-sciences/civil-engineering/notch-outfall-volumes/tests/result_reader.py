"""graded_pair(text) -> {"volume_m3": float, "peak_15min_L_s": float}; only these two values are ever converted.

Two passes. The first replaces every JSON string by a numbered placeholder token (quotes and backslash escapes
handled by a regular expression), so that the second pass sees brackets, commas and colons only where they are
structural. The second pass cuts the outer object at depth-1 commas; each piece must be <placeholder> : <raw value>.
Placeholders are decoded back into names; a graded name must appear once, with a raw value that is a bare JSON
number and a finite float. Inner objects, arrays and every other member are never parsed.
"""
import json
import math
import re

GRADED = ("volume_m3", "peak_15min_L_s")
STRING = re.compile(r'"(?:[^"\\]|\\.)*"', re.S)
NUMBER = re.compile(r"-?(?:0|[1-9][0-9]*)(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?")
TOKEN = re.compile(r"\x00(\d+)\x00")


class BadResult(ValueError):
    pass


def _mask(text):
    strings = []

    def keep(match):
        strings.append(match.group(0))
        return "\x00%d\x00" % (len(strings) - 1)

    masked = STRING.sub(keep, text)
    if '"' in masked or "\x00" in STRING.sub("", text):
        raise BadResult("unterminated string")
    return masked, strings


def _top_level(masked):
    body = masked.strip()
    if not body.startswith("{"):
        raise BadResult("the result must be one JSON object")
    depth, start, pieces = 0, 1, []
    for i, ch in enumerate(body):
        if ch in "{[":
            depth += 1
        elif ch in "}]":
            depth -= 1
            if depth < 0:
                raise BadResult("unbalanced brackets")
            if depth == 0:
                if body[i + 1:].strip():
                    raise BadResult("text after the object")
                pieces.append(body[start:i])
                return [p for p in pieces if p.strip()]
        elif ch == "," and depth == 1:
            pieces.append(body[start:i])
            start = i + 1
    raise BadResult("object not closed")


def graded_pair(text):
    masked, strings = _mask(text)
    out = {}
    for piece in _top_level(masked):
        head, sep, raw = piece.partition(":")
        name_token = TOKEN.fullmatch(head.strip())
        if not sep or not name_token:
            raise BadResult("member without a string name")
        name = json.loads(strings[int(name_token.group(1))])
        if name not in GRADED:
            continue
        if name in out:
            raise BadResult(f"{name} twice")
        raw = raw.strip()
        if not NUMBER.fullmatch(raw):
            raise BadResult(f"{name} is not a bare number")
        try:
            value = float(raw)
        except (OverflowError, ValueError):
            raise BadResult(f"{name} cannot be converted")
        if not math.isfinite(value):
            raise BadResult(f"{name} is not finite")
        out[name] = value
    missing = [g for g in GRADED if g not in out]
    if missing:
        raise BadResult("missing " + ", ".join(missing))
    return out
