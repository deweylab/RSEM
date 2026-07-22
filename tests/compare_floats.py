#!/usr/bin/env python3
"""Compare two whitespace-delimited numeric files with relative-error tolerance.

Used instead of exact `diff` for outputs (e.g. .theta) that carry tiny
floating-point differences across platforms even on a correct run. Relative
error only, so near-zero values may still differ on a large relative swing.
"""
import sys
import math


def rows(path):
    with open(path) as f:
        return [line.split() for line in f]


def close(a, b, rtol=1e-4):
    try:
        return math.isclose(float(a), float(b), rel_tol=rtol)
    except ValueError:
        return a == b  # non-numeric tokens must match exactly


def main():
    if len(sys.argv) != 3:
        sys.exit(f"usage: {sys.argv[0]} <file_a> <file_b>")
    a_path, b_path = sys.argv[1], sys.argv[2]
    a_rows, b_rows = rows(a_path), rows(b_path)
    if len(a_rows) != len(b_rows):
        sys.exit(f"{a_path} vs {b_path}: line count {len(a_rows)} != {len(b_rows)}")
    for n, (a, b) in enumerate(zip(a_rows, b_rows), 1):
        if len(a) != len(b):
            sys.exit(f"{a_path}:{n} vs {b_path}:{n}: field count {len(a)} != {len(b)}")
        for c, (va, vb) in enumerate(zip(a, b), 1):
            if not close(va, vb):
                sys.exit(f"{a_path}:{n}:{c}: {va} != {vb} (outside tolerance)")
    print(f"OK (within tolerance): {a_path} vs {b_path}")


if __name__ == "__main__":
    main()
