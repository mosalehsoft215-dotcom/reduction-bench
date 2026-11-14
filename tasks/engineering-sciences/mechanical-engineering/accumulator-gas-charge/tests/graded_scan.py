"""graded_pair(text) -> {name: float} for the two graded members of a result object; nothing else is decoded.

class Scanner keeps a position in the text and a bracket count.  At bracket count 1 it reads a member name (decoded
only to compare it), a colon and a value; a bare scalar value is kept as raw text, any string or container value is
stepped over character by character (string escapes respected, brackets counted, never decoded).  After the closing
brace only blanks may follow.  Graded names must occur once, as finite JSON numbers.
"""
import json
import math
import re

GRADED = ("n2_precharge_bar", "gas_volume_litres")
NUMBER = re.compile(r"-?(?:0|[1-9][0-9]*)(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?")
BLANK = " \t\r\n"
STOPS = BLANK + ",]}:[{\""


class BadResult(ValueError):
    pass


class Scanner:
    def __init__(self, text):
        self.text, self.at = text, 0

    def peek(self):
        while self.at < len(self.text) and self.text[self.at] in BLANK:
            self.at += 1
        return self.text[self.at] if self.at < len(self.text) else ""

    def take(self, wanted):
        if self.peek() != wanted:
            raise BadResult(f"expected {wanted!r} at {self.at}")
        self.at += 1

    def string_end(self):
        """index just past the string that starts at self.at"""
        i = self.at + 1
        while i < len(self.text):
            if self.text[i] == "\\":
                i += 2
            elif self.text[i] == '"':
                return i + 1
            else:
                i += 1
        raise BadResult("string not closed")

    def skip_container(self):
        depth = 0
        while self.at < len(self.text):
            ch = self.text[self.at]
            if ch == '"':
                self.at = self.string_end()
                continue
            if ch in "[{":
                depth += 1
            elif ch in "]}":
                depth -= 1
                if depth == 0:
                    self.at += 1
                    return
            self.at += 1
        raise BadResult("container not closed")

    def scalar(self):
        start = self.at
        while self.at < len(self.text) and self.text[self.at] not in STOPS:
            self.at += 1
        if self.at == start:
            raise BadResult(f"value expected at {start}")
        return self.text[start:self.at]

    def members(self):
        self.take("{")
        if self.peek() == "}":
            self.at += 1
            return []
        found = []
        while True:
            if self.peek() != '"':
                raise BadResult("member name expected")
            end = self.string_end()
            try:
                name = json.loads(self.text[self.at:end])
            except ValueError:
                raise BadResult("bad member name")
            self.at = end
            self.take(":")
            head = self.peek()
            if head and head in "[{":
                self.skip_container()
                found.append((name, None))
            elif head == '"':
                self.at = self.string_end()
                found.append((name, None))
            else:
                found.append((name, self.scalar()))
            if self.peek() == ",":
                self.at += 1
                continue
            self.take("}")
            if self.peek():
                raise BadResult("text after the object")
            return found


def graded_pair(text):
    raw = {}
    for name, tok in Scanner(text).members():
        if name in GRADED:
            if name in raw:
                raise BadResult(name + " given twice")
            raw[name] = tok
    result = {}
    for name in GRADED:
        tok = raw.get(name)
        if tok is None or not NUMBER.fullmatch(tok):
            raise BadResult(name + " missing or not a number")
        try:
            value = float(tok)
        except (OverflowError, ValueError):
            raise BadResult(name + " not finite")
        if not math.isfinite(value):
            raise BadResult(name + " not finite")
        result[name] = value
    return result
