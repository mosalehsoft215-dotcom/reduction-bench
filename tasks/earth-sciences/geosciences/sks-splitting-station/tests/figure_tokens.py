"""Finds the top-level figures of an answer from its text, without parsing anything else.

lexemes(text) yields the JSON lexemes one by one. figures(text, names) walks them with a depth counter and a list of
open brackets and returns {name: [(lexeme kind, lexeme)...]} for the wanted names, or None if the text is not one valid
JSON object. A number of any length, any nesting, NaN, Infinity or a repeated name under another key changes nothing.
"""
import json
import re

PIECE = re.compile(r'\s*(?:(?P<s>"(?:[^"\\\x00-\x1f]|\\["\\/bfnrt]|\\u[0-9a-fA-F]{4})*")'
                   r'|(?P<n>-?(?:0|[1-9]\d*)(?:\.\d+)?(?:[eE][-+]?\d+)?)'
                   r'|(?P<c>true|false|null|NaN|-?Infinity)|(?P<p>[\[\]{},:]))')


def lexemes(text):
    end = len(text.rstrip())
    place = 0
    while place < end:
        got = PIECE.match(text, place)
        if got is None:
            yield "bad", ""
            return
        place = got.end()
        yield got.lastgroup, got.group(got.lastgroup)


def figures(text, names):
    out = {name: [] for name in names}
    opened = []
    expecting = {"value"}
    current = None
    begun = False
    for kind, word in lexemes(text):
        if kind == "bad":
            return None
        if not begun:
            if word != "{":
                return None
            begun = True
            opened.append("}")
            expecting = {"key", "close"}
            continue
        if not opened:
            return None
        if "close" in expecting and word == opened[-1]:
            opened.pop()
            expecting = {"comma", "close"}
        elif "comma" in expecting and word == ",":
            expecting = {"key"} if opened[-1] == "}" else {"value"}
        elif "key" in expecting and kind == "s":
            if len(opened) == 1:
                current = json.loads(word)
            expecting = {"colon"}
        elif "colon" in expecting and word == ":":
            expecting = {"value"}
        elif "value" in expecting and word in ("{", "["):
            if len(opened) == 1 and current in out:
                out[current].append(("container", word))
            opened.append("}" if word == "{" else "]")
            expecting = {"key", "close"} if word == "{" else {"value", "close"}
        elif "value" in expecting and kind in ("s", "n", "c"):
            if len(opened) == 1 and current in out:
                out[current].append((kind, word))
            expecting = {"comma", "close"}
        else:
            return None
    return out if begun and not opened else None
