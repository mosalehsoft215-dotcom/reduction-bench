"""force_values(text) -> {key: float}: the two graded members of a result object, nothing else decoded.

Method: one pass marks every character that lies inside a JSON string (escapes honoured).  The nesting depth of
each remaining character is the running count of opening minus closing brackets (numpy cumulative sums).  Commas
and colons at depth 1 cut the outer object into members; a member's name is decoded only to compare it, and a value
is taken as raw text only when its name is graded.  Graded names must occur once each, as finite JSON numbers.
"""
import json
import math
import re

import numpy as np

KEYS = ("axial_before_kn", "axial_after_kn")
NUMBER = re.compile(r"-?(?:0|[1-9][0-9]*)(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?")


class Unreadable(ValueError):
    pass


def string_mask(text):
    mask = np.zeros(len(text), dtype=bool)
    inside = escaped = False
    for i, ch in enumerate(text):
        if inside:
            mask[i] = True
            if escaped:
                escaped = False
            elif ch == "\\":
                escaped = True
            elif ch == '"':
                inside = False
        elif ch == '"':
            inside = True
            mask[i] = True
    if inside:
        raise Unreadable("string not closed")
    return mask


def outer_members(text):
    body = text.strip()
    if not body.startswith("{") or not body.endswith("}"):
        raise Unreadable("not one JSON object")
    chars = np.frombuffer(body.encode("utf-32-le"), dtype=np.uint32)
    quoted = string_mask(body)
    step = np.zeros(len(body), dtype=np.int64)
    for ch, delta in (("{", 1), ("[", 1), ("}", -1), ("]", -1)):
        step[(chars == ord(ch)) & ~quoted] += delta
    depth_after = np.cumsum(step)
    if depth_after.min() < 0 or depth_after[-1] != 0 or np.any(depth_after[:-1] == 0):
        raise Unreadable("brackets do not balance into one object")
    opens = [i for i, c in enumerate(body) if c in "[{" and not quoted[i]]
    closes = [i for i, c in enumerate(body) if c in "]}" and not quoted[i]]
    if len(opens) != len(closes):
        raise Unreadable("brackets do not match")
    stack = []
    for i in sorted(opens + closes):
        if body[i] in "[{":
            stack.append(body[i])
        elif not stack or {"]": "[", "}": "{"}[body[i]] != stack.pop():
            raise Unreadable("brackets do not match")
    level_one = (depth_after == 1) & ~quoted
    cuts = [i for i in range(1, len(body) - 1) if level_one[i] and body[i] == ","]
    inner = body[1:-1].strip()
    if not inner:
        return []
    found = []
    for a, b in zip([0] + cuts, cuts + [len(body) - 1]):
        piece_start = a + 1
        colon = next((j for j in range(piece_start, b) if level_one[j] and body[j] == ":"), None)
        if colon is None:
            raise Unreadable("member without a colon")
        raw_name = body[piece_start:colon].strip()
        if not (raw_name.startswith('"') and raw_name.endswith('"')) or len(raw_name) < 2:
            raise Unreadable("member name is not a string")
        try:
            name = json.loads(raw_name)
        except ValueError:
            raise Unreadable("bad member name")
        value = body[colon + 1:b].strip()
        if not value:
            raise Unreadable("member without a value")
        found.append((name, value))
    return found


def force_values(text):
    raw = {}
    for name, value in outer_members(text):
        if name in KEYS:
            if name in raw:
                raise Unreadable(name + " appears twice")
            raw[name] = value
    values = {}
    for name in KEYS:
        token = raw.get(name)
        if token is None or not NUMBER.fullmatch(token):
            raise Unreadable(name + " absent or not a plain number")
        try:
            number = float(token)
        except (OverflowError, ValueError):
            raise Unreadable(name + " is not finite")
        if not math.isfinite(number):
            raise Unreadable(name + " is not finite")
        values[name] = number
    return values
