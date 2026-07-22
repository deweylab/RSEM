#!/usr/bin/env python3
"""Compare two whitespace-delimited numeric files.

Used instead of exact `diff` for outputs that differ across platforms only in
ways that aren't real regressions:
  default    relative-error tolerance, for values with tiny float drift (.theta)
  --exact    exact match per token, for values that should be identical (.model)

Both modes treat nan and -nan as equal (the NaN sign bit is compiler-specific).
"""
import sys
import math


def rows(path):
    with open(path) as f:
        return [line.split() for line in f]


def equal(a, b, exact, rtol=1e-4):
    try:
        fa, fb = float(a), float(b)
    except ValueError:
        return a == b  # non-numeric tokens must match exactly
    if math.isnan(fa) and math.isnan(fb):
        return True  # nan vs -nan: same value, sign bit differs across compilers
    if exact:
        return fa == fb
    return math.isclose(fa, fb, rel_tol=rtol)


def main():
    args = [a for a in sys.argv[1:] if a != "--exact"]
    exact = "--exact" in sys.argv[1:]
    if len(args) != 2:
        sys.exit(f"usage: {sys.argv[0]} [--exact] <file_a> <file_b>")
    a_path, b_path = args
    a_rows, b_rows = rows(a_path), rows(b_path)
    if len(a_rows) != len(b_rows):
        sys.exit(f"{a_path} vs {b_path}: line count {len(a_rows)} != {len(b_rows)}")
    for n, (a, b) in enumerate(zip(a_rows, b_rows), 1):
        if len(a) != len(b):
            sys.exit(f"{a_path}:{n} vs {b_path}:{n}: field count {len(a)} != {len(b)}")
        for c, (va, vb) in enumerate(zip(a, b), 1):
            if not equal(va, vb, exact):
                sys.exit(f"{a_path}:{n}:{c}: {va} != {vb}")
    print(f"OK: {a_path} vs {b_path}")


if __name__ == "__main__":
    main()
