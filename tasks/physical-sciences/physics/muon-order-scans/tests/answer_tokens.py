"""Reads the two values of a submitted result without materialising the JSON document.

A single pass over the text produces structural events (object and array brackets, member names, scalars), keeping
the open brackets on a plain Python list so that very deep nesting cannot exhaust the interpreter stack. The
consumer validates the grammar and keeps `local_field_mt` and `neel_temperature_k` only where they are members of
the top-level object: each has to be present exactly once and be a JSON number. All other content is validated
and discarded; NaN and the infinities, which Python's json writes, are tolerated anywhere except in those two.
"""
import json
import math
import re

WANTED = ("local_field_mt", "neel_temperature_k")
NUM = re.compile(r"-?(?:0|[1-9]\d*)(?:\.\d+)?(?:[eE][+-]?\d+)?")
STR = re.compile(r'"(?:[^"\\\x00-\x1f]|\\["\\/bfnrt]|\\u[0-9a-fA-F]{4})*"')
WORDS = ("true", "false", "null", "NaN", "Infinity", "-Infinity")
BLANK = " \t\r\n"


class Rejected(ValueError):
    pass


def events(text):
    i, n = 0, len(text)

    def skip(i):
        while i < n and text[i] in BLANK:
            i += 1
        return i

    def scalar(i):
        m = STR.match(text, i)
        if m:
            return ("string", m.group()), m.end()
        for w in WORDS:
            if text.startswith(w, i):
                return ("word", w), i + len(w)
        m = NUM.match(text, i)
        if m:
            return ("number", m.group()), m.end()
        raise Rejected("no value at %d" % i)

    stack = []
    i = skip(i)
    if i >= n or text[i] != "{":
        raise Rejected("the outer value is not an object")
    want_value = True
    while True:
        i = skip(i)
        if i >= n:
            raise Rejected("text ends early")
        ch = text[i]
        if want_value:
            if ch in "{[":
                stack.append(ch)
                yield ("open", ch, len(stack))
                i += 1
                i = skip(i)
                if i < n and text[i] == ("}" if ch == "{" else "]"):
                    stack.pop()
                    yield ("close", ch, len(stack) + 1)
                    i += 1
                    want_value = False
                    if not stack:
                        break
                    continue
                if ch == "{":
                    m = STR.match(text, i)
                    if not m:
                        raise Rejected("member name expected")
                    yield ("name", m.group(), len(stack))
                    i = skip(m.end())
                    if i >= n or text[i] != ":":
                        raise Rejected("colon expected")
                    i += 1
                continue
            item, i = scalar(i)
            yield ("scalar", item, len(stack))
            want_value = False
            continue
        top = stack[-1]
        if ch == ",":
            i = skip(i + 1)
            if top == "{":
                m = STR.match(text, i)
                if not m:
                    raise Rejected("member name expected")
                yield ("name", m.group(), len(stack))
                i = skip(m.end())
                if i >= n or text[i] != ":":
                    raise Rejected("colon expected")
                i += 1
            want_value = True
            continue
        if ch == ("}" if top == "{" else "]"):
            stack.pop()
            yield ("close", top, len(stack) + 1)
            i += 1
            if not stack:
                break
            continue
        raise Rejected("comma or closing bracket expected at %d" % i)
    if skip(i) != n:
        raise Rejected("text after the outer object")


def graded(text):
    found, pending, seen = {}, None, set()
    for kind, item, depth in events(text):
        if depth != 1:
            continue
        if kind == "name":
            name = json.loads(item)
            pending = name if name in WANTED else None
            if pending is not None:
                if pending in seen:
                    raise Rejected(pending + " appears twice")
                seen.add(pending)
        elif pending is not None:
            if kind != "scalar" or item[0] != "number":
                raise Rejected(pending + " is not a number")
            found[pending] = item[1]
            pending = None
        elif kind == "scalar" or kind == "open":
            pending = None
    for k in WANTED:
        if k not in found:
            raise Rejected("missing " + k)
        if not math.isfinite(float(found[k])):
            raise Rejected(k + " is not finite")
    return found


def read(path):
    try:
        with open(path, "rb") as fh:
            return graded(fh.read().decode("utf-8"))
    except (OSError, UnicodeDecodeError, ValueError):
        return None
