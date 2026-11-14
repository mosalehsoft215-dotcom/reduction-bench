"""read_graded(text): the two graded numbers of an output document, nothing else decoded.

A pushdown walk over the characters: a list holds the open brackets, a flag marks string context, and at stack
height 1 (inside the outer object) the walk alternates between member names and values.  Only member names are
decoded; a value is copied as raw text when its name is graded and is otherwise passed over.  Graded names must
appear once at the top level with a finite JSON number.
"""
import json
import math
import re
import sys

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

GRADED = ("precharge_bar_20c", "gas_volume_l")
NUMERIC = re.compile(r"-?(?:0|[1-9]\d*)(?:\.\d+)?(?:[eE][+-]?\d+)?\Z")


class Rejected(ValueError):
    pass


def _top_level_members(text):
    stack, in_str, esc = [], False, False
    name_start = value_start = None
    expect = "name"
    members = []
    name = None
    for pos, ch in enumerate(text):
        if in_str:
            if esc:
                esc = False
            elif ch == "\\":
                esc = True
            elif ch == '"':
                in_str = False
                if len(stack) == 1 and expect == "name_end":
                    name = json.loads(text[name_start:pos + 1])
                    expect = "colon"
            continue
        if ch.isspace():
            continue
        depth = len(stack)
        if ch == '"':
            in_str = True
            if depth == 1 and expect == "name":
                name_start, expect = pos, "name_end"
            elif depth == 1 and expect == "value":
                value_start, expect = pos, "in_value"
            continue
        if depth == 0:
            if ch != "{" or members or name is not None:
                raise Rejected("not one JSON object")
            stack.append("{")
            continue
        if depth == 1 and expect == "colon":
            if ch != ":":
                raise Rejected("colon expected")
            expect = "value"
            continue
        if depth == 1 and expect == "value":
            value_start, expect = pos, "in_value"
        if ch in "[{":
            stack.append(ch)
        elif ch in "]}":
            if not stack or {"]": "[", "}": "{"}[ch] != stack[-1]:
                raise Rejected("brackets do not match")
            if len(stack) == 1:
                if expect == "in_value":
                    members.append((name, text[value_start:pos].strip()))
                elif expect != "name" or name is not None:
                    raise Rejected("object ends early")
                stack.pop()
                expect = "done"
                continue
            stack.pop()
        elif ch == "," and len(stack) == 1:
            if expect != "in_value":
                raise Rejected("misplaced comma")
            members.append((name, text[value_start:pos].strip()))
            name, expect = None, "name"
        elif expect == "done":
            raise Rejected("text after the object")
    if stack or in_str or expect != "done":
        raise Rejected("object not closed")
    return members


def read_graded(text):
    raw = {}
    for name, value in _top_level_members(text):
        if name in GRADED:
            if name in raw:
                raise Rejected(f"{name} repeated")
            raw[name] = value
    out = {}
    for name in GRADED:
        value = raw.get(name)
        if value is None or not NUMERIC.match(value):
            raise Rejected(f"{name} absent or not a number")
        try:
            number = float(value)
        except (OverflowError, ValueError):
            raise Rejected(f"{name} not finite")
        if not math.isfinite(number):
            raise Rejected(f"{name} not finite")
        out[name] = number
    return out
