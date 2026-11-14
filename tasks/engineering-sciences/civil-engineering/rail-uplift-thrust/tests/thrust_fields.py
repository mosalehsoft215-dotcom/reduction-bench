"""pick_graded(text) -> {key: float}: the two graded members of a result object, read without decoding the rest.

A table-driven automaton: every character is mapped to a class (quote, backslash, opener, closer, colon, comma,
blank, other) and STEP[(state, class)] gives the next state plus a mark.  A nesting counter is kept for values that
open a container, so nested text is only counted, never decoded.  Marks note where a member name or a bare scalar
starts and ends; names are decoded alone to compare them.  Graded keys must occur once, as finite JSON numbers.
"""
import json
import math
import re

KEYS = ("destress_temp_c", "rail_thrust_kn")
JSON_NUMBER = re.compile(r"-?(?:0|[1-9]\d*)(?:\.\d+)?(?:[eE][+-]?\d+)?")
CLASS = {'"': "q", "\\": "b", "{": "o", "[": "o", "}": "c", "]": "c", ":": "k", ",": "m",
         " ": "s", "\t": "s", "\r": "s", "\n": "s"}

# state -> {class: (next state, mark)}; a missing entry is an error.  "*" matches any class.
STEP = {
    "start": {"s": ("start", None), "o": ("key?", "open")},
    "key?": {"s": ("key?", None), "q": ("key", "key_from"), "c": ("done", "close_empty")},
    "key": {"b": ("key\\", None), "q": ("colon", "key_to"), "*": ("key", None)},
    "key\\": {"*": ("key", None)},
    "colon": {"s": ("colon", None), "k": ("value?", None)},
    "value?": {"s": ("value?", None), "q": ("text", None), "o": ("deep", "dive"), "b": None, "c": None, "k": None,
               "m": None, "*": ("bare", "bare_from")},
    "bare": {"s": ("next", "bare_to"), "m": ("key!", "bare_to"), "c": ("done", "bare_to_close"), "*": ("bare", None)},
    "text": {"b": ("text\\", None), "q": ("next", "skipped"), "*": ("text", None)},
    "text\\": {"*": ("text", None)},
    "deep": {"q": ("deep_text", None), "o": ("deep", "dive"), "c": ("deep", "rise"), "*": ("deep", None)},
    "deep_text": {"b": ("deep_text\\", None), "q": ("deep", None), "*": ("deep_text", None)},
    "deep_text\\": {"*": ("deep_text", None)},
    "next": {"s": ("next", None), "m": ("key!", None), "c": ("done", "close")},
    "key!": {"s": ("key!", None), "q": ("key", "key_from")},
    "done": {"s": ("done", None)},
}


class Unreadable(ValueError):
    pass


def members(text):
    state, depth, found = "start", 0, []
    name_from = value_from = 0
    name = None
    for pos, ch in enumerate(text):
        table = STEP[state]
        cls = CLASS.get(ch, "x")
        move = table.get(cls, table.get("*"))
        if move is None:
            raise Unreadable(f"unexpected {ch!r} at {pos}")
        state, mark = move
        if mark == "key_from":
            name_from = pos
        elif mark == "key_to":
            try:
                name = json.loads(text[name_from:pos + 1])
            except ValueError:
                raise Unreadable("bad name")
        elif mark == "bare_from":
            value_from = pos
        elif mark in ("bare_to", "bare_to_close"):
            found.append((name, text[value_from:pos]))
        elif mark == "skipped":
            found.append((name, None))
        elif mark == "dive":
            depth += 1
        elif mark == "rise":
            depth -= 1
            if depth == 0:
                found.append((name, None))
                state = "next"
    if state != "done":
        raise Unreadable("object not complete")
    return found


def pick_graded(text):
    seen = {}
    for name, raw in members(text):
        if name in KEYS:
            if name in seen:
                raise Unreadable(f"{name} repeated")
            seen[name] = raw
    picked = {}
    for key in KEYS:
        raw = seen.get(key)
        if raw is None or JSON_NUMBER.fullmatch(raw) is None:
            raise Unreadable(f"{key} absent or not a number")
        try:
            value = float(raw)
        except (OverflowError, ValueError):
            raise Unreadable(f"{key} not finite")
        if not math.isfinite(value):
            raise Unreadable(f"{key} not finite")
        picked[key] = value
    return picked
