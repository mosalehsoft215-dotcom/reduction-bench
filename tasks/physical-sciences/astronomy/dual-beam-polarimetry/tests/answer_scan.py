"""Pulls the two graded numbers out of a result file without building any JSON value.

A regular-expression tokenizer cuts the text into JSON tokens; a small loop then checks the grammar with a stack
of open brackets and, while the stack holds only the outer object, notes the members `target_p_percent` and
`target_angle_deg`. Every other member is grammar-checked and dropped, so depth and size cost only time. NaN,
Infinity and -Infinity are accepted as values (Python's json module writes them), but never for a graded member.

    graded(text) -> {"target_p_percent": "<number text>", "target_angle_deg": "<number text>"}  or raises Bad
"""
import json
import math
import re

KEYS = ("target_p_percent", "target_angle_deg")
TOKEN = re.compile(r'''
    (?P<ws>[ \t\n\r]+)
  | (?P<str>"(?:[^"\\\x00-\x1f]|\\["\\/bfnrt]|\\u[0-9a-fA-F]{4})*")
  | (?P<num>-?(?:0|[1-9][0-9]*)(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)
  | (?P<lit>true|false|null|NaN|-Infinity|Infinity)
  | (?P<punct>[{}\[\],:])
''', re.X)


class Bad(ValueError):
    pass


def tokens(text):
    pos, end = 0, len(text)
    while pos < end:
        m = TOKEN.match(text, pos)
        if not m:
            raise Bad("unreadable text at %d" % pos)
        pos = m.end()
        if m.lastgroup != "ws":
            yield m.lastgroup, m.group()


def graded(text):
    stream = tokens(text)
    first = next(stream, None)
    if first != ("punct", "{"):
        raise Bad("the outer value is not an object")
    opened = ["{"]
    state = "first_key"            # first_key | key | colon | value | first_value | after
    current = None                 # graded key whose value comes next (outer object only)
    found = {}
    for kind, tok in stream:
        if not opened:
            raise Bad("text after the outer object")
        top = opened[-1]
        if state in ("first_key", "key"):
            if kind == "punct" and tok == "}" and state == "first_key":
                opened.pop(); state = "after"; continue
            if kind != "str":
                raise Bad("a member name was expected")
            if len(opened) == 1:
                name = json.loads(tok)
                current = name if name in KEYS else None
                if current in found:
                    raise Bad(current + " occurs twice")
            state = "colon"
        elif state == "colon":
            if (kind, tok) != ("punct", ":"):
                raise Bad("a colon was expected")
            state = "value"
        elif state in ("value", "first_value"):
            if kind == "punct" and tok == "]" and state == "first_value":
                opened.pop(); state = "after"; continue
            outer = len(opened) == 1 and top == "{"
            if kind == "punct" and tok in "{[":
                if outer and current:
                    raise Bad(current + " is not a number")
                opened.append(tok)
                state = "first_key" if tok == "{" else "first_value"
                continue
            if kind == "punct":
                raise Bad("a value was expected")
            if outer and current:
                if kind != "num":
                    raise Bad(current + " is not a number")
                found[current] = tok
            if outer:
                current = None
            state = "after"
        else:                       # after a value or a closed container
            closer = "}" if top == "{" else "]"
            if (kind, tok) == ("punct", ","):
                state = "key" if top == "{" else "value"
            elif (kind, tok) == ("punct", closer):
                opened.pop()
                state = "after"
            else:
                raise Bad("a comma or a closing bracket was expected")
    if opened:
        raise Bad("the text ends inside a value")
    for k in KEYS:
        if k not in found:
            raise Bad("missing " + k)
        if not math.isfinite(float(found[k])):
            raise Bad(k + " is not finite")
    return found


def from_file(path):
    try:
        with open(path, "rb") as fh:
            return graded(fh.read().decode("utf-8"))
    except (OSError, UnicodeDecodeError, ValueError):
        return None
