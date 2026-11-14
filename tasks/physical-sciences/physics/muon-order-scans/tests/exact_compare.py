"""Exact recount of a window check with Python fractions, no floating point: the reported number's text and the
truth's repr are turned into exact fractions through decimal.Decimal, and |reported - truth| <= width is decided on
those. Exponents beyond 400 in magnitude are out of any physical range and count as outside."""
from decimal import Decimal, InvalidOperation
from fractions import Fraction


def exact(text):
    d = Decimal(text.strip())
    if not d.is_finite() or abs(d.as_tuple().exponent) > 400:
        raise ValueError(text)
    return Fraction(d)


def inside(reported, truth, width):
    try:
        return abs(exact(reported) - exact(truth)) <= exact(width)
    except (InvalidOperation, ValueError):
        return False
