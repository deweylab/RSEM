#!/usr/bin/env python3
"""Ablation check for option-coverage test cases.

For each named case, run rsem-calculate-expression once with the case's full
argument list, then once more per flag with that flag reverted to RSEM's own
documented default. If a reverted run produces byte-identical output to the
full run, that flag isn't actually exercising anything on this dataset, and
the test case is not real coverage for it.

This is a discussion/starting-point script, not wired into `make test` yet --
see the NOTE below on fragment-length-mean/sd before trusting it as-is.
"""
import filecmp
import subprocess
import sys
from pathlib import Path

TEST_READS_PAIRED = ["--paired-end", "tests/data/reads_1.fastq", "tests/data/reads_2.fastq"]
TEST_READS_SINGLE = ["tests/data/reads_1.fastq"]
REFERENCE = "tests/gold/paired_end/bowtie/reference/my_ref"
REFERENCE_SE = "tests/gold/single_end/bowtie/reference/my_ref"

# name -> {base: full arg list (minus reads/reference), ablate: {flag: default_or_None}}
# "default_or_None": the RSEM-documented default value to substitute; None means the flag
# is boolean and should just be dropped entirely for the ablated run.
CASES = {
    "bowtie_custom": {
        "base": ["--bowtie-n", "3", "--bowtie-e", "200", "--bowtie-m", "5", "--seed-length", "28"],
        "reads": TEST_READS_PAIRED,
        "reference": REFERENCE,
        "ablate": {
            "--bowtie-n": "2",        # RSEM/Bowtie default
            "--bowtie-e": "99999999",  # RSEM/Bowtie default
            "--bowtie-m": "200",       # RSEM/Bowtie default
            "--seed-length": "25",     # RSEM default
        },
    },
    "bowtie2_custom": {
        "base": ["--bowtie2", "--bowtie2-mismatch-rate", "0.05", "--bowtie2-k", "5", "--bowtie2-sensitivity-level", "very_fast"],
        "reads": TEST_READS_PAIRED,
        "reference": "tests/gold/paired_end/bowtie2/reference/my_ref",
        "ablate": {
            "--bowtie2-mismatch-rate": "0.1",       # RSEM default
            "--bowtie2-k": "200",                    # RSEM default
            "--bowtie2-sensitivity-level": "sensitive",  # RSEM default
        },
    },
    "gibbs_sampling": {
        "base": ["--seed", "0", "--single-cell-prior", "--calc-pme", "--calc-ci",
                 "--gibbs-burnin", "20", "--gibbs-number-of-samples", "800",
                 "--gibbs-sampling-gap", "2", "--ci-credibility-level", "0.80",
                 "--ci-number-of-samples-per-count-vector", "30"],
        "reads": TEST_READS_PAIRED,
        "reference": REFERENCE,
        "ablate": {
            "--gibbs-burnin": "200",                         # RSEM default
            "--gibbs-number-of-samples": "1000",             # RSEM default
            "--gibbs-sampling-gap": "1",                     # RSEM default
            "--ci-credibility-level": "0.95",                # RSEM default
            "--ci-number-of-samples-per-count-vector": "50", # RSEM default
            "--single-cell-prior": None,                      # boolean -> drop it
        },
    },
    "fragment_modeling": {
        "base": ["--fragment-length-mean", "400", "--fragment-length-sd", "50",
                 "--fragment-length-min", "100", "--fragment-length-max", "700",
                 "--estimate-rspd", "--num-rspd-bins", "40"],
        "reads": TEST_READS_SINGLE,
        "reference": REFERENCE_SE,
        "ablate": {
            "--fragment-length-min": "1",     # RSEM default
            "--fragment-length-max": "1000",  # RSEM default
            "--num-rspd-bins": "20",          # RSEM default
            # NOTE: --fragment-length-mean's RSEM default (-1) *disables* fragment-length
            # modeling entirely, and --fragment-length-sd's default depends on --mean being
            # set. Reverting either to its literal default changes the shape of the run, not
            # just the value -- needs a decision (e.g. compare against a mean/sd pair that's
            # merely "different" rather than "off") before this pair can ablate meaningfully.
        },
    },
}


def run(outdir: Path, args, reads, reference, sample="my_sample"):
    outdir.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["./rsem-calculate-expression", "-p", "1", *args, *reads, reference, str(outdir / sample)],
        check=True,
    )


def ablate(args, flag, default):
    out = list(args)
    i = out.index(flag)
    if default is None:
        del out[i]
    else:
        out[i + 1] = default
    return out


def main():
    failures = []
    for case, spec in CASES.items():
        base_dir = Path(f"tests/output/coverage/{case}/base")
        run(base_dir, spec["base"], spec["reads"], spec["reference"])
        for flag, default in spec["ablate"].items():
            variant_dir = Path(f"tests/output/coverage/{case}/without_{flag.strip('-')}")
            run(variant_dir, ablate(spec["base"], flag, default), spec["reads"], spec["reference"])
            if filecmp.cmp(base_dir / "my_sample.isoforms.results", variant_dir / "my_sample.isoforms.results", shallow=False):
                failures.append(f"{case}: {flag} produced IDENTICAL output when reverted to default -- not exercised")
    if failures:
        sys.exit("Option coverage check FAILED:\n" + "\n".join(f" - {f}" for f in failures))
    print("Every ablated flag changed output -- coverage confirmed.")


if __name__ == "__main__":
    main()
