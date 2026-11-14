"""Reads `stress_free_temp_c` and `axial_force_kn` from a result document; nothing else in it is decoded.

Two passes.  The first pass masks every JSON string (a regular expression finds them, escapes included) and then
records the nesting depth at every character, so commas and colons that matter are exactly those at depth 1.  The
second pass cuts the outer object at its depth-1 commas into members, splits each member at its first depth-1 colon,
decodes only the member name, and keeps the raw value text of the two wanted names.  A wanted name must occur once
and hold a JSON number that is finite as a binary float; a malformed outer object is refused.
"""
import json
import math
import re
import sys

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

WANTED = ("stress_free_temp_c", "axial_force_kn")
STRING = re.compile(r'"(?:[^"\\]|\\.)*"', re.S)
NUMBER = re.compile(r"-?(?:0|[1-9]\d*)(?:\.\d+)?(?:[eE][+-]?\d+)?")


class Unusable(ValueError):
    pass


def _masked(text):
    chars = list(text)
    for m in STRING.finditer(text):
        for k in range(m.start() + 1, m.end() - 1):
            chars[k] = "s"
    masked = "".join(chars)
    if masked.count('"') % 2:
        raise Unusable("unbalanced quotes")
    return masked


def _members(text):
    body = text.strip()
    if not (body.startswith("{") and body.endswith("}")):
        raise Unusable("not one JSON object")
    inner = body[1:-1]
    masked = _masked(inner)
    pieces, depth, start = [], 0, 0
    for k, c in enumerate(masked):
        if c in "[{":
            depth += 1
        elif c in "]}":
            depth -= 1
            if depth < 0:
                raise Unusable("unbalanced brackets")
        elif c == "," and depth == 0:
            pieces.append((inner[start:k], masked[start:k]))
            start = k + 1
    if depth != 0:
        raise Unusable("unbalanced brackets")
    if inner[start:].strip():
        pieces.append((inner[start:], masked[start:]))
    elif pieces:
        raise Unusable("trailing comma")
    for raw, mask in pieces:
        cut = mask.find(":")
        if cut < 0:
            raise Unusable("member without a colon")
        name_text = raw[:cut].strip()
        if not (name_text.startswith('"') and name_text.endswith('"')):
            raise Unusable("member name is not a string")
        value = raw[cut + 1:].strip()
        if not value:
            raise Unusable("member without a value")
        yield json.loads(name_text), value


def pick(text):
    seen = {}
    for name, value in _members(text):
        if name in WANTED:
            if name in seen:
                raise Unusable(f"{name} appears twice")
            seen[name] = value
    result = {}
    for name in WANTED:
        raw = seen.get(name)
        if raw is None or not NUMBER.fullmatch(raw):
            raise Unusable(f"{name} missing or not a number")
        try:
            value = float(raw)
        except (OverflowError, ValueError):
            raise Unusable(f"{name} not finite")
        if not math.isfinite(value):
            raise Unusable(f"{name} not finite")
        result[name] = value
    return result
