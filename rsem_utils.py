#!/usr/bin/env python3
"""
RSEM utility functions - Python port of rsem_perl_utils.pm
"""

import os
import sys
import subprocess

VERSION = "RSEM v1.3.1"  # Update version info here
SAMTOOLS = "samtools-1.3"  # If update to another version of SAMtools, need to change this

ALLELE_TITLE = [
    "allele_id", "transcript_id", "gene_id", "length", "effective_length",
    "expected_count", "TPM", "FPKM", "AlleleIsoPct", "AlleleGenePct",
    "posterior_mean_count", "posterior_standard_deviation_of_count",
    "pme_TPM", "pme_FPKM", "AlleleIsoPct_from_pme_TPM", "AlleleGenePct_from_pme_TPM",
    "TPM_ci_lower_bound", "TPM_ci_upper_bound", "TPM_coefficient_of_quartile_variation",
    "FPKM_ci_lower_bound", "FPKM_ci_upper_bound", "FPKM_coefficient_of_quartile_variation",
]

TRANSCRIPT_TITLE = [
    "transcript_id", "gene_id", "length", "effective_length", "expected_count",
    "TPM", "FPKM", "IsoPct", "posterior_mean_count", "posterior_standard_deviation_of_count",
    "pme_TPM", "pme_FPKM", "IsoPct_from_pme_TPM",
    "TPM_ci_lower_bound", "TPM_ci_upper_bound", "TPM_coefficient_of_quartile_variation",
    "FPKM_ci_lower_bound", "FPKM_ci_upper_bound", "FPKM_coefficient_of_quartile_variation",
]

GENE_TITLE = [
    "gene_id", "transcript_id(s)", "length", "effective_length", "expected_count",
    "TPM", "FPKM", "posterior_mean_count", "posterior_standard_deviation_of_count",
    "pme_TPM", "pme_FPKM",
    "TPM_ci_lower_bound", "TPM_ci_upper_bound", "TPM_coefficient_of_quartile_variation",
    "FPKM_ci_lower_bound", "FPKM_ci_upper_bound", "FPKM_coefficient_of_quartile_variation",
]


def run_command(cmd, err_msg=None):
    """Execute a shell command. Exit on failure."""
    print(cmd)
    try:
        result = subprocess.run(cmd, shell=True)
        if result.returncode != 0:
            err = err_msg + "\n" if err_msg else ""
            err += f'"{cmd}" failed! Please check if you provide correct parameters/options for the pipeline!\n'
            print(err, file=sys.stderr)
            sys.exit(-1)
    except OSError as e:
        first_word = cmd.split()[0] if cmd.split() else cmd
        print(f"{first_word} : {e}", file=sys.stderr)
        print('Please check if you have compiled the associated codes by typing related "make" commands and/or made related executables ready to use.', file=sys.stderr)
        sys.exit(-1)
    print()


def collect_results(result_type, inp_f, out_f):
    """
    Transpose and add headers to RSEM result files.
    result_type: "allele", "isoform", or "gene"
    """
    try:
        with open(inp_f) as f:
            results = [line.rstrip("\n").split("\t") for line in f]
    except OSError:
        print(f"Fail to open file {inp_f}!", file=sys.stderr)
        sys.exit(-1)

    try:
        with open(out_f, "w") as f:
            n = len(results)
            m = len(results[0]) if results else 0

            if result_type == "allele":
                titles = ALLELE_TITLE
            elif result_type == "isoform":
                titles = TRANSCRIPT_TITLE
            elif result_type == "gene":
                titles = GENE_TITLE
            else:
                print("A bug on 'collectResults' is detected!", file=sys.stderr)
                sys.exit(-1)

            f.write("\t".join(titles[i] for i in range(n)) + "\n")
            for i in range(m):
                row = [results[j][i] for j in range(n)]
                f.write("\t".join(row) + "\n")
    except OSError:
        print(f"Fail to create file {out_f}!", file=sys.stderr)
        sys.exit(-1)


def show_version_info():
    """Print version and exit."""
    print(f"Current version: {VERSION}")
    sys.exit(0)


def get_samtools():
    """Return SAMtools directory name."""
    return SAMTOOLS


def has_poly_a(filename):
    """Check if file indicates polyA (fullLen < totLen)."""
    with open(filename) as f:
        line = f.readline().rstrip("\n")
    full_len, tot_len = line.split()
    return int(full_len) < int(tot_len)
