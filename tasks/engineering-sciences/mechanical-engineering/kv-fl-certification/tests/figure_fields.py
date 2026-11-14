"""pick(text) -> {"kv_rated_m3h": float, "fl_rated": float}; nothing else in the text is ever converted.

A generator walks the text once and yields (depth, piece) for every outer-level member: it tracks string state with
escapes and a stack of open brackets, and emits the raw text of a member each time a comma or the final brace closes
it at depth 1. pick() then takes the name up to the first colon outside strings, decodes only that name, and for the
two graded names requires a bare JSON number that converts to a finite float, each name once. Exactly one object,
nothing after it.
"""
import json
import math
import re

GRADED = ("kv_rated_m3h", "fl_rated")
BARE_NUMBER = re.compile(r"-?(?:0|[1-9][0-9]*)(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?")
CLOSE_OF = {"{": "}", "[": "]"}


class Unreadable(ValueError):
    pass


def outer_members(text):
    s = text.strip()
    if not s.startswith("{"):
        raise Unreadable("figures must be one JSON object")
    stack, in_string, escape, mark = [], False, False, 1
    for i, c in enumerate(s):
        if in_string:
            escape, in_string = (False, True) if escape else (c == "\\", c != '"')
            continue
        if c == '"':
            in_string = True
        elif c in CLOSE_OF:
            stack.append(CLOSE_OF[c])
        elif c in "]}":
            if not stack or stack.pop() != c:
                raise Unreadable("unbalanced brackets")
            if not stack:
                if s[i + 1:].strip():
                    raise Unreadable("text after the object")
                yield s[mark:i]
                return
        elif c == "," and len(stack) == 1:
            yield s[mark:i]
            mark = i + 1
    raise Unreadable("object never closed")


def name_and_value(member):
    in_string, escape, depth = False, False, 0
    for i, c in enumerate(member):
        if in_string:
            escape, in_string = (False, True) if escape else (c == "\\", c != '"')
            continue
        if c == '"':
            in_string = True
        elif c in "{[":
            depth += 1
        elif c in "]}":
            depth -= 1
        elif c == ":" and depth == 0:
            return member[:i].strip(), member[i + 1:].strip()
    raise Unreadable("member without a value")


def pick(text):
    seen = {}
    members = [m for m in outer_members(text)]
    if len(members) == 1 and not members[0].strip():
        members = []
    for member in members:
        raw_name, raw_value = name_and_value(member)
        if len(raw_name) < 2 or raw_name[0] != '"' or raw_name[-1] != '"':
            raise Unreadable("member name is not a string")
        try:
            name = json.loads(raw_name)
        except ValueError:
            raise Unreadable("member name unreadable")
        if name not in GRADED:
            continue
        if name in seen:
            raise Unreadable(f"{name} appears twice")
        if not BARE_NUMBER.fullmatch(raw_value):
            raise Unreadable(f"{name} is not a bare number")
        try:
            number = float(raw_value)
        except (OverflowError, ValueError):
            raise Unreadable(f"{name} does not convert")
        if not math.isfinite(number):
            raise Unreadable(f"{name} is not finite")
        seen[name] = number
    absent = [g for g in GRADED if g not in seen]
    if absent:
        raise Unreadable("absent: " + ", ".join(absent))
    return seen
