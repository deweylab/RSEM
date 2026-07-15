#!/usr/bin/env python3
"""Compare two whitespace-delimited numeric files, allowing small relative-error tolerance.

Used in place of exact `diff` for RSEM outputs (e.g. .theta files) that can carry
tiny floating-point differences across platforms/compilers even on a correct run.

Relative error only (no absolute floor): matches the convergence criterion RSEM's own
EM algorithm uses (see the paper's Methods -- "stopped when all theta_i with value >=
1e-7 have a relative change of less than 1e-3"). Values very close to zero can therefore
still trip this check on a large relative swing between two tiny numbers; that's a known
tradeoff of dropping the absolute-tolerance floor, not a bug.
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
