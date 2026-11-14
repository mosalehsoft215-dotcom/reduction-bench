"""Outer members of a JSON answer, taken from its text without building any value.

Walker.feed() pushes one lexeme at a time through a small table of states, with a list as the nesting stack, so a number
of any length, nesting of any depth, NaN, Infinity or a repeated name in a key that is not graded changes nothing.
outer_members(text) gives {key: [(sort, lexeme), ...]} for one valid JSON object, or None for anything else.
"""
import json
import re

LEXEME = re.compile(r'[ \t\r\n]*("(?:[^"\\\x00-\x1f]|\\["\\/bfnrt]|\\u[0-9a-fA-F]{4})*"'
                    r'|-?(?:0|[1-9][0-9]*)(?:\.[0-9]+)?(?:[eE][-+]?[0-9]+)?'
                    r'|true|false|null|NaN|-?Infinity|[\[\]{},:])')


def sort_of(lexeme):
    if lexeme[0] == '"':
        return "string"
    if lexeme[0].isdigit() or (lexeme[0] == "-" and lexeme != "-Infinity"):
        return "number"
    return lexeme if lexeme in "{}[],:" else "constant"


class Walker:
    def __init__(self):
        self.open = []
        self.state = "start"
        self.name = None
        self.seen = {}

    def note(self, sort, lexeme):
        if len(self.open) == 1:
            self.seen.setdefault(self.name, []).append((sort, lexeme))

    def enter(self, lexeme):
        self.note("container", lexeme)
        self.open.append(lexeme)
        self.state = "first key" if lexeme == "{" else "first item"

    def leave(self):
        self.open.pop()
        self.state = "after"

    def feed(self, lexeme):
        sort = sort_of(lexeme)
        state = self.state
        if state == "start":
            if lexeme != "{":
                return False
            self.open.append(lexeme)
            self.state = "first key"
        elif state in ("first key", "key"):
            if lexeme == "}" and state == "first key":
                self.leave()
            elif sort == "string":
                if len(self.open) == 1:
                    self.name = json.loads(lexeme)
                self.state = "colon"
            else:
                return False
        elif state == "colon":
            if lexeme != ":":
                return False
            self.state = "item"
        elif state in ("first item", "item"):
            if lexeme == "]" and state == "first item":
                self.leave()
            elif lexeme in ("{", "["):
                self.enter(lexeme)
            elif sort in ("string", "number", "constant"):
                self.note(sort, lexeme)
                self.state = "after"
            else:
                return False
        else:
            if not self.open:
                return False
            closer = "}" if self.open[-1] == "{" else "]"
            if lexeme == closer:
                self.leave()
            elif lexeme == ",":
                self.state = "key" if self.open[-1] == "{" else "item"
            else:
                return False
        return True


def outer_members(text):
    walker, place, stop = Walker(), 0, len(text.rstrip(" \t\r\n"))
    while place < stop:
        found = LEXEME.match(text, place)
        if found is None or not walker.feed(found.group(1)):
            return None
        place = found.end()
    return walker.seen if walker.state == "after" and not walker.open else None
