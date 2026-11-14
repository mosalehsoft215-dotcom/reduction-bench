"""Reads a result as text and keeps only the members of its outer object; both judges go through this file.

outer_pairs(text) gives [(key, kind, raw token)] per member, or None when the text is not one JSON object. The pass uses
an explicit list for nesting and never builds a value, so whatever other keys hold (numbers of any length, any nesting,
NaN, Infinity, repeated names) has no effect on the keys that are judged.
"""
import json
import re

TOKEN = re.compile(r'\s*(?:(?P<text>"(?:[^"\\\x00-\x1f]|\\(?:["\\/bfnrt]|u[0-9a-fA-F]{4}))*")'
                   r'|(?P<number>-?(?:0|[1-9]\d*)(?:\.\d+)?(?:[eE][-+]?\d+)?)'
                   r'|(?P<word>true|false|null|NaN|-?Infinity)|(?P<mark>[{}\[\],:]))')


def outer_pairs(text):
    pairs_seen, nest, expect, key, at = [], [], "value", None, 0
    end = len(text.rstrip())
    while at < end:
        hit = TOKEN.match(text, at)
        if hit is None:
            return None
        at, kind, raw = hit.end(), hit.lastgroup, hit.group(hit.lastgroup)
        if expect in ("value", "value or close"):
            if raw == "]" and expect == "value or close":
                nest.pop()
                expect = "next"
            elif kind != "mark":
                if len(nest) == 1:
                    pairs_seen.append((key, kind, raw))
                expect = "next"
            elif raw in "{[":
                if len(nest) == 1:
                    pairs_seen.append((key, "container", raw))
                nest.append(raw)
                expect = "key or close" if raw == "{" else "value or close"
            else:
                return None
        elif expect in ("key", "key or close"):
            if raw == "}" and expect == "key or close":
                nest.pop()
                expect = "next"
            elif kind == "text":
                if len(nest) == 1:
                    key = json.loads(raw)
                expect = "colon"
            else:
                return None
        elif expect == "colon":
            if raw != ":":
                return None
            expect = "value"
        else:
            if not nest:
                return None
            if raw == ",":
                expect = "key" if nest[-1] == "{" else "value"
            elif raw == ("}" if nest[-1] == "{" else "]"):
                nest.pop()
            else:
                return None
    fine = not nest and expect == "next" and text.lstrip().startswith("{")
    return pairs_seen if fine else None
