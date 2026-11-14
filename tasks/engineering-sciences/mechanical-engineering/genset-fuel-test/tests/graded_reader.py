"""Pulls `sfc75_g_per_kwh` and `sfc100_g_per_kwh` out of a result without decoding the whole JSON value.

The text is cut into tokens by one regular expression and walked by a small grammar machine that keeps only a list of
open containers, so arbitrarily deep nesting is cheap. The two names count only as members of the outermost object;
each must be present exactly once with a JSON number as its value. Every other member is checked for grammar and then
forgotten, and NaN or the infinities are tolerated there because Python's json module writes them.
"""
import math
import re

NAMES = ("sfc75_g_per_kwh", "sfc100_g_per_kwh")
TOKEN = re.compile(r'\s*(?:(?P<s>"(?:[^"\\\x00-\x1f]|\\["\\/bfnrt]|\\u[0-9a-fA-F]{4})*")'
                   r'|(?P<n>-?(?:0|[1-9][0-9]*)(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)(?![0-9A-Za-z.])'
                   r'|(?P<w>true|false|null|NaN|-?Infinity)(?![0-9A-Za-z])|(?P<p>[{}\[\]:,]))')


class Malformed(Exception):
    pass


class Walker:
    def __init__(self, text):
        self.text, self.pos = text, 0
        self.found = {}
        self._obj_open = False
        self._array_open = False

    def tokens(self):
        while True:
            m = TOKEN.match(self.text, self.pos)
            if not m:
                if self.text[self.pos:].strip():
                    raise Malformed(f"unexpected text at {self.pos}")
                return
            self.pos = m.end()
            kind = m.lastgroup
            yield kind, m.group(kind)

    def run(self):
        stack = []                       # 'o' or 'a' per open container
        expect = "value"                 # value | key | colon | after
        pending = None                   # graded name waiting for its value at depth 1
        seen_root = False
        for kind, tok in self.tokens():
            if expect == "key":
                if kind == "s":
                    self._obj_open = False
                    pending = tok[1:-1] if (len(stack) == 1 and tok[1:-1] in NAMES) else None
                    if pending and pending in self.found:
                        raise Malformed(f"{pending} given twice")
                    expect = "colon"
                    continue
                if tok == "}" and stack[-1] == "o" and self._obj_open:
                    stack.pop()
                    self._obj_open = False
                    expect = "after" if stack else "end"
                    continue
                raise Malformed("member name expected")
            if expect == "colon":
                if tok != ":":
                    raise Malformed("':' expected")
                expect = "value"
                continue
            if expect == "value":
                if seen_root and not stack:
                    raise Malformed("text after the root value")
                seen_root = True
                if pending is not None:
                    if kind != "n":
                        raise Malformed(f"{pending} is not a number")
                    self.found[pending] = tok
                    pending = None
                if tok in ("{", "["):
                    stack.append("o" if tok == "{" else "a")
                    self._obj_open = tok == "{"
                    expect = "key" if tok == "{" else "value_or_close"
                    if expect == "value_or_close":
                        expect = "value"
                        self._array_open = True
                    continue
                if tok == "]" and stack and stack[-1] == "a" and getattr(self, "_array_open", False):
                    stack.pop()
                    self._array_open = False
                    expect = "after"
                    continue
                if kind in ("s", "n", "w"):
                    self._array_open = False
                    expect = "after"
                    continue
                raise Malformed("value expected")
            if expect == "after":
                if not stack:
                    raise Malformed("text after the root value")
                if tok == ",":
                    expect = "key" if stack[-1] == "o" else "value"
                    self._array_open = False
                    continue
                if (tok == "}" and stack[-1] == "o") or (tok == "]" and stack[-1] == "a"):
                    stack.pop()
                    expect = "after" if stack else "end"
                    continue
                raise Malformed("',' or closing bracket expected")
            if expect == "end":
                raise Malformed("text after the root value")
        if stack or not seen_root:
            raise Malformed("unterminated JSON text")
        return self.found


def values(text):
    found = Walker(text).run()
    if set(found) != set(NAMES):
        raise Malformed("missing " + ", ".join(sorted(set(NAMES) - set(found))))
    for name, tok in found.items():
        if not math.isfinite(float(tok)):
            raise Malformed(f"{name} is not finite")
    return found


def load(path):
    """The two values as their number text, or None when the result is unusable."""
    try:
        return values(path.read_bytes().decode("utf-8"))
    except (OSError, UnicodeDecodeError, Malformed, ValueError):
        return None
