"""Pulls `ic77_a` and `n77` out of an answer document without decoding anything else in it.

The text is cut into JSON tokens by one regular expression (strings with their escapes, the six structural
characters, and bare atoms).  Only tokens at nesting depth 1 matter: a member name, its colon, and - when the name
is graded - the single atom that follows.  Everything deeper is counted past, so members that are not graded may hold
anything of any size or depth.  The outer object must be well formed, each graded name must occur once at depth 1,
and its value must be a JSON number that stays finite as a binary float.
"""
import json
import math
import re
import sys

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

WANTED = ("ic77_a", "n77")
TOKENS = re.compile(r'\s*(?:("(?:[^"\\]|\\.)*")|([{}\[\]:,])|([^\s{}\[\]:,"]+))', re.S)
JSON_NUMBER = re.compile(r"-?(?:0|[1-9]\d*)(?:\.\d+)?(?:[eE][+-]?\d+)?")


class Refused(ValueError):
    pass


def tokens(text):
    pos = 0
    while True:
        m = TOKENS.match(text, pos)
        if not m or m.end() == pos:
            if text[pos:].strip():
                raise Refused("unreadable text")
            return
        pos = m.end()
        yield m.group(1), m.group(2), m.group(3)


def graded_members(text):
    depth = 0
    state = "start"            # start -> name -> colon -> value -> sep -> name ... -> done
    current = None
    raw = {}
    for string, mark, atom in tokens(text):
        if state == "done":
            raise Refused("text after the object")
        if depth == 0:
            if mark != "{":
                raise Refused("not a JSON object")
            depth, state = 1, "name"
            continue
        if depth > 1:                                      # inside a member value that is not looked at
            if mark in ("[", "{"):
                depth += 1
            elif mark in ("]", "}"):
                depth -= 1
                if depth == 1:
                    state = "sep"
            continue
        if state == "name":
            if mark == "}" and not raw and current is None:
                state, depth = "done", 0
                continue
            if string is None:
                raise Refused("member name expected")
            current, state = json.loads(string), "colon"
        elif state == "colon":
            if mark != ":":
                raise Refused("colon expected")
            state = "value"
        elif state == "value":
            if mark in ("[", "{"):
                if current in WANTED:
                    raise Refused(f"{current} is not a number")
                depth, state = 2, "nested"
                continue
            if mark is not None:
                raise Refused("value expected")
            if current in WANTED:
                if current in raw:
                    raise Refused(f"{current} given twice")
                raw[current] = atom if atom is not None else string
            state = "sep"
        elif state == "sep":
            if mark == ",":
                state = "name"
            elif mark == "}":
                state, depth = "done", 0
            else:
                raise Refused("comma or closing brace expected")
    if state != "done":
        raise Refused("object not closed")
    values = {}
    for name in WANTED:
        token = raw.get(name)
        if token is None:
            raise Refused(f"{name} missing")
        if not JSON_NUMBER.fullmatch(token):
            raise Refused(f"{name} is not a JSON number")
        try:
            number = float(token)
        except (OverflowError, ValueError):
            raise Refused(f"{name} is not finite")
        if not math.isfinite(number):
            raise Refused(f"{name} is not finite")
        values[name] = number
    return values
