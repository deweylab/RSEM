#!/usr/bin/env python

from sys import argv, exit

allele_title = ["allele_id", "transcript_id", "gene_id", "length", "effective_length", "expected_count", "TPM", "FPKM", "AlleleIsoPct", "AlleleGenePct", "posterior_mean_count", "posterior_standard_deviation_of_count", "pme_TPM", "pme_FPKM", "AlleleIsoPct_from_pme_TPM", "AlleleGenePct_from_pme_TPM", "TPM_ci_lower_bound", "TPM_ci_upper_bound", "TPM_coefficient_of_quartile_variation", "FPKM_ci_lower_bound", "FPKM_ci_upper_bound", "FPKM_coefficient_of_quartile_variation"]

transcript_title = ["transcript_id", "gene_id", "length", "effective_length", "expected_count", "TPM", "FPKM", "IsoPct", "posterior_mean_count", "posterior_standard_deviation_of_count", "pme_TPM", "pme_FPKM", "IsoPct_from_pme_TPM", "TPM_ci_lower_bound", "TPM_ci_upper_bound", "TPM_coefficient_of_quartile_variation", "FPKM_ci_lower_bound", "FPKM_ci_upper_bound", "FPKM_coefficient_of_quartile_variation"]

gene_title = ["gene_id", "transcript_id(s)", "length", "effective_length", "expected_count", "TPM", "FPKM", "posterior_mean_count", "posterior_standard_deviation_of_count", "pme_TPM", "pme_FPKM", "TPM_ci_lower_bound", "TPM_ci_upper_bound", "TPM_coefficient_of_quartile_variation", "FPKM_ci_lower_bound", "FPKM_ci_upper_bound", "FPKM_coefficient_of_quartile_variation"]

if len(argv) != 4:
    print("python extract_results.py type input output")
    exit(-1)

with open(argv[2]) as fin:
    values = []
    for line in fin:
        values.append(line.strip().split('\t'))

title = None
if argv[1] == "gene":
    title = gene_title
elif argv[1] == "isoform":
    title = transcript_title
else:
    title = allele_title

nc = len(values)
M = len(values[0])

with open(argv[3], "w") as fout:
    fout.write("\t".join(title[:nc]) + "\n")
    for i in xrange(M):
        for j in xrange(nc - 1):
            fout.write(values[j][i] + "\t")
        fout.write(values[nc - 1][i] + "\n")
