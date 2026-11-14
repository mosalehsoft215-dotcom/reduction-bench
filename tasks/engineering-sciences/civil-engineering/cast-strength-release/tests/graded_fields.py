"""pick(text) -> the two graded numbers of a target object; nothing else is decoded.

A transition table drives the scan: each character is classed (blank, quote, opener, closer, colon, comma, other) and
TABLE[(mode, class)] names the step to take.  Strings and nested containers met as values are jumped over by the
`hop` step (escapes respected, brackets counted, nothing decoded).  Member names are decoded only to be compared; bare
values are cut out as raw text.  Graded names must appear once at the top level as finite JSON numbers.
"""
import json, math, re

GRADED = ("strength_at_release_mpa", "time_to_10mpa_h")
NUMBER = re.compile(r"-?(?:0|[1-9]\d*)(?:\.\d+)?(?:[eE][+-]?\d+)?")
CLASS = {" ": "blank", "\t": "blank", "\r": "blank", "\n": "blank", '"': "quote", "{": "open", "[": "open",
         "}": "close", "]": "close", ":": "colon", ",": "comma"}


class Unreadable(ValueError):
    pass


def hop(text, at):
    """index just past the string or container starting at `at`"""
    depth, i = 0, at
    while i < len(text):
        ch = text[i]
        if ch == '"':
            i += 1
            while i < len(text) and text[i] != '"':
                i += 2 if text[i] == "\\" else 1
            if i >= len(text):
                raise Unreadable("unterminated string")
            if depth == 0:
                return i + 1
        elif ch in "[{":
            depth += 1
        elif ch in "]}":
            depth -= 1
            if depth == 0:
                return i + 1
        i += 1
    raise Unreadable("unterminated value")


def scan(text):
    s = {"mode": "begin", "name": None, "start": 0, "pairs": []}

    def opened(i):
        if text[i] != "{":
            raise Unreadable("not a JSON object")
        s["mode"] = "first"
        return i + 1

    def name(i):
        end = hop(text, i)
        try:
            s["name"] = json.loads(text[i:end])
        except ValueError:
            raise Unreadable("bad member name")
        s["mode"] = "colon"
        return end

    def colon(i):
        s["mode"] = "value"
        return i + 1

    def skipped(i):
        s["pairs"].append((s["name"], None))
        s["mode"] = "next"
        return hop(text, i)

    def bare(i):
        s["start"], s["mode"] = i, "bare"
        return i + 1

    def bare_more(i):
        return i + 1

    def bare_end(i):
        s["pairs"].append((s["name"], text[s["start"]:i]))
        s["mode"] = "next"
        return i

    def comma(i):
        s["mode"] = "member"
        return i + 1

    def closed(i):
        if text[i] != "}":
            raise Unreadable("brackets do not match")
        s["mode"] = "done"
        return i + 1

    def blank(i):
        return i + 1

    table = {("begin", "open"): opened, ("first", "quote"): name, ("member", "quote"): name,
             ("first", "close"): closed, ("colon", "colon"): colon, ("value", "quote"): skipped,
             ("value", "open"): skipped, ("value", "other"): bare, ("bare", "other"): bare_more,
             ("next", "comma"): comma, ("next", "close"): closed}
    for mode in ("bare",):
        for cls in ("blank", "quote", "open", "close", "colon", "comma"):
            table[(mode, cls)] = bare_end
    for mode in ("begin", "first", "member", "colon", "value", "next", "done"):
        table[(mode, "blank")] = blank
    i = 0
    while i < len(text):
        cls = CLASS.get(text[i], "other")
        step = table.get((s["mode"], cls))
        if step is None:
            raise Unreadable(f"unexpected {text[i]!r} while in {s['mode']}")
        i = step(i)
    if s["mode"] != "done":
        raise Unreadable("object not closed")
    return s["pairs"]


def pick(text):
    raw = {}
    for name, value in scan(text):
        if name in GRADED:
            if name in raw:
                raise Unreadable(name + " appears twice")
            raw[name] = value
    out = {}
    for name in GRADED:
        value = raw.get(name)
        if value is None or not NUMBER.fullmatch(value):
            raise Unreadable(name + " missing or not a JSON number")
        try:
            number = float(value)
        except (OverflowError, ValueError):
            raise Unreadable(name + " not finite")
        if not math.isfinite(number):
            raise Unreadable(name + " not finite")
        out[name] = number
    return out
