"""Validate a JSON object iteratively and retain its requested numeric tokens."""
import re
from json.decoder import scanstring

SPACE = re.compile(r"[ \t\r\n]*")
SCALAR = re.compile(r"(?:-?(?:0|[1-9][0-9]*)(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?|true|false|null)")


def figures(text, names):
    results = {name: [] for name in names}
    at = SPACE.match(text).end()
    if text[at:at + 1] != "{":
        return None
    at += 1
    # Each frame carries its closer and one of three expectations:
    # initial item or close, required item, separator or close.
    stack = [("}", "initial")]
    try:
        while stack:
            closer, phase = stack[-1]
            at = SPACE.match(text, at).end()
            mark = text[at]
            if mark == closer and phase != "required":
                stack.pop()
                at += 1
                continue
            if phase == "separator":
                if mark != ",":
                    return None
                stack[-1] = (closer, "required")
                at += 1
                continue
            key = None
            if closer == "}":
                if mark != '"':
                    return None
                key, at = scanstring(text, at + 1)
                at = SPACE.match(text, at).end()
                if text[at] != ":":
                    return None
                at = SPACE.match(text, at + 1).end()
                mark = text[at]
            slot = results.get(key) if len(stack) == 1 else None
            stack[-1] = (closer, "separator")
            if mark == "{" or mark == "[":
                stack.append(("}" if mark == "{" else "]", "initial"))
                at += 1
                item = ("container", mark)
            elif mark == '"':
                _, at = scanstring(text, at + 1)
                item = ("string", "")
            else:
                token = SCALAR.match(text, at)
                if token is None:
                    return None
                word = token.group()
                at = token.end()
                item = ("literal" if word in ("true", "false", "null") else "n", word)
            if slot is not None:
                slot.append(item)
    except (IndexError, ValueError):
        return None
    return results if SPACE.fullmatch(text, at) else None
