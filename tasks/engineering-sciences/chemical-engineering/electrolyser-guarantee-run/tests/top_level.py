"""Reads the members of the outermost JSON object of an answer without decoding the values nested inside it.

Nested objects and arrays are stepped over by counting brackets outside strings, so no depth of nesting can exhaust
a recursion limit. Strings are matched by the JSON string grammar; scalars must be JSON numbers or one of true,
false, null, NaN, Infinity, -Infinity. Only the outermost member names matter to the judges; a graded name must
occur exactly once there and its value must be written as a JSON number.
"""
import re

STRING = re.compile(r'"(?:[^"\\\x00-\x1f]|\\["\\/bfnrt]|\\u[0-9a-fA-F]{4})*"')
NUMBER = re.compile(r"-?(?:0|[1-9][0-9]*)(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?")
WORD = re.compile(r"true|false|null|NaN|-?Infinity")
SPACE = " \t\r\n"


class BadAnswer(Exception):
    pass


def _space(text, i):
    while i < len(text) and text[i] in SPACE:
        i += 1
    return i


def _value_end(text, i):
    if i >= len(text):
        raise BadAnswer("value missing")
    if text[i] == '"':
        m = STRING.match(text, i)
        if not m:
            raise BadAnswer("bad string")
        return m.end()
    if text[i] in "[{":
        depth = 0
        while i < len(text):
            ch = text[i]
            if ch == '"':
                m = STRING.match(text, i)
                if not m:
                    raise BadAnswer("bad string")
                i = m.end()
                continue
            if ch in "[{":
                depth += 1
            elif ch in "]}":
                depth -= 1
                if depth == 0:
                    return i + 1
            i += 1
        raise BadAnswer("unterminated container")
    m = NUMBER.match(text, i) or WORD.match(text, i)
    if not m:
        raise BadAnswer("bad scalar")
    return m.end()


def outer_members(text):
    i = _space(text, 0)
    if i >= len(text) or text[i] != "{":
        raise BadAnswer("not an object")
    i = _space(text, i + 1)
    members = []
    if i < len(text) and text[i] == "}":
        i += 1
    else:
        while True:
            m = STRING.match(text, i)
            if not m:
                raise BadAnswer("member name expected")
            name = m.group()[1:-1]
            i = _space(text, m.end())
            if i >= len(text) or text[i] != ":":
                raise BadAnswer("':' expected")
            start = _space(text, i + 1)
            end = _value_end(text, start)
            members.append((name, text[start:end]))
            i = _space(text, end)
            if i < len(text) and text[i] == ",":
                i = _space(text, i + 1)
                continue
            if i < len(text) and text[i] == "}":
                i += 1
                break
            raise BadAnswer("',' or '}' expected")
    if _space(text, i) != len(text):
        raise BadAnswer("text after the object")
    return members


def graded_raw(text, names):
    """{name: number text} for the graded names; BadAnswer when one is missing, repeated or not a number."""
    found = {}
    for name, raw in outer_members(text):
        if name in names:
            if name in found:
                raise BadAnswer(f"{name} repeated")
            if not NUMBER.fullmatch(raw):
                raise BadAnswer(f"{name} is not a JSON number")
            found[name] = raw
    missing = [n for n in names if n not in found]
    if missing:
        raise BadAnswer("missing " + ", ".join(missing))
    return found
