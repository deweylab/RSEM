#!/usr/bin/env python

from sys import argv, exit
from math import floor
import pysam

try:
	xrange
except NameError:
	xrange = range

if len(argv) != 4:
	print("Usage: python correctbackground.py reference_name sample_name number_of_bins")
	exit(-1)


reference_fasta = argv[1] + ".idx.fa"
sampleToken = argv[2][argv[2].rfind('/') + 1 : ]
model_file = "{0}.stat/{1}.model".format(argv[2], sampleToken)
expression_file = "{0}.isoforms.results".format(argv[2])
bam_file = "{0}.transcript.bam".format(argv[2])
gene_expression_file = "{0}.genes.results".format(argv[2])
out_expression_file = "{0}_corrected.isoforms.results".format(argv[2])
out_gene_expression_file = "{0}_corrected.genes.results".format(argv[2])


def find_insert_size_mode(model_file):
	with open(model_file) as fin:
		for i in range(4):
			next(fin)
		lb = int(next(fin).strip().split()[0])
		values = [float(x) for x in next(fin).strip().split()]
		return lb + 1 + [i for i in range(len(values)) if values[i] == max(values)][0]


isize = find_insert_size_mode(model_file)
n = int(argv[3])
background = [0.0] * n
foreground = [0.0] * n
GCfactor = [0.0] * n
n_GC = 0 

def load_ref(reference_fasta):
	refs = []
	names = []
	with open(reference_fasta) as fin:
		for line in fin:
			names.append(line[1:-1])
			refs.append(next(fin).strip())

	print("load_ref is done.")

	return names, refs

def isGC(c):
	return 1 if (c == 'G' or c == 'C') else 0

def get_category(refseq, start, length):
	global n_GC
	n_GC = sum([1 if isGC(x) else 0 for x in refseq[start : (start + length)]])
	category = int((n_GC * 1.0 / length) * n)
	return category if category < n else n - 1

def get_cat_nGC(nGC):
	category = int((n_GC * 1.0 / isize) * n)
	return category if category < n else n - 1



def calc_background(names, refs, expression_file):
	global background, n_GC

	tid = 0
	with open(expression_file) as fin:
		next(fin)
		for line in fin:
			fields = line.split()
			tpm = float(fields[5])
			assert names[tid] == fields[0]
			tlen = len(refs[tid])
			background[get_category(refs[tid], 0, min(tlen, isize))] += tpm
			for pos in xrange(isize, tlen):
				n_GC = n_GC + isGC(refs[tid][pos]) - isGC(refs[tid][pos - isize])
				background[get_cat_nGC(n_GC)] += tpm

			tid += 1
			if tid % 5000 == 0:
				print("FIN {0}.".format(tid))

	denom = sum(background)
	background = [x / denom for x in background]

	print("calc_background is done.")



def buildMap(references, names):
	assert len(references) <= len(names)

	local_map = {}
	for i in xrange(len(names)):
		local_map[names[i]] = i

	ref_id_to_tid = []
	for ref_name in references:
		ref_name = ref_name.split()[0]
		tid = local_map.get(ref_name, -1)
		assert tid >= 0
		ref_id_to_tid.append(tid)

	return ref_id_to_tid

def calc_foreground(names, refs, bam_file):
	global foreground

	sam_in = pysam.AlignmentFile(bam_file, "rb")

	ref_id_to_tid = buildMap(sam_in.references, names)

	cnt = 0
	for aread in sam_in:
		assert aread.is_paired
		if not aread.is_unmapped and not aread.is_reverse:
			frac = aread.get_tag("ZW")
			if frac > 0.0:
				category = get_category(refs[ref_id_to_tid[aread.reference_id]], aread.reference_start, abs(aread.template_length))
				foreground[category] += frac

		cnt += 1
		if cnt % 1000000 == 0:
			print("FIN {0}.".format(cnt))

	sam_in.close()

	denom = sum(foreground)
	foreground = [x / denom for x in foreground]

	print("calc_foreground is done.")



def calcGCfactor(foreground, background):
	return [x / y if y > 0.0 else 0.0 for x, y in zip(foreground, background)]

def calcEffLen(sequence):
	global n_GC

	length = len(sequence)
	efflen = GCfactor[get_category(sequence, 0, min(isize, length))]
	for pos in xrange(isize, length):
		n_GC = n_GC + isGC(sequence[pos]) - isGC(sequence[pos - isize])
		efflen += GCfactor[get_cat_nGC(n_GC)]

	return efflen

def generate_expression_files(refs, expression_file, gene_expression_file, out_expression_file, out_gene_expression_file):
	efflens = []
	exprs = []

	isopct = []
	gefflens = []
	gexprs = []

	cgid = ""
	cgexpr = 0.0
	cisopct = []

	tid = 0
	with open(expression_file) as fin:
		next(fin)
		for line in fin:
			fields = line.strip().split()
			efflen = calcEffLen(refs[tid])
			expr = float(fields[4]) / efflen if efflen > 0.0 else 0.0

			if cgid != fields[1]:
				if cgid != "":
					cisopct = [x / cgexpr if cgexpr > 0.0 else 0.0 for x in cisopct]
					isopct.extend(cisopct)

					gexprs.append(cgexpr)

					nt = len(cisopct)
					gefflens.append(sum([x * y for x, y in zip(cisopct, efflens[-nt : ])]))

				cgid = fields[1]
				cgexpr = 0.0
				cisopct = []

			efflens.append(efflen)
			exprs.append(expr)

			cgexpr += expr
			cisopct.append(expr)

			tid += 1
			if tid % 5000 == 0:
				print("FIN {0}.".format(tid))


	if cgid != "":
		cisopct = [x / cgexpr if cgexpr > 0.0 else 0.0 for x in cisopct]
		isopct.extend(cisopct)

		gexprs.append(cgexpr)

		nt = len(cisopct)
		gefflens.append(sum([x * y for x, y in zip(cisopct, efflens[-nt : ])]))

	denom = sum(exprs)
	exprs = [x / denom * 1e6 for x in exprs]
	factor = 1e9 / sum([x * y for x, y in zip(exprs, efflens)])
	fpkms = [x * factor for x in exprs]

	denom = sum(gexprs)
	gexprs = [x / denom * 1e6 for x in gexprs]
	gfpkms = [x * factor for x in gexprs]

	with open(expression_file) as fin, open(out_expression_file, "w") as fout:
		 fields = next(fin).strip().split()
		 fout.write("\t".join(fields[:8]) + "\n")

		 tid = 0
		 for line in fin:
		 	fields = line.strip().split()
		 	fields[3] = "{0:.2f}".format(efflens[tid])
		 	fields[5] = "{0:.2f}".format(exprs[tid])
		 	fields[6] = "{0:.2f}".format(fpkms[tid])
		 	fields[7] = "{0:.2f}".format(isopct[tid] * 1e2)
		 	fout.write("\t".join(fields[:8]) + "\n")

		 	tid += 1

	with open(gene_expression_file) as fin, open(out_gene_expression_file, "w") as fout:
		fields = next(fin).strip().split()
		fout.write("\t".join(fields[:7]) + "\n")

		gid = 0
		for line in fin:
			fields = line.strip().split()
			fields[3] = "{0:.2f}".format(gefflens[gid])
			fields[5] = "{0:.2f}".format(gexprs[gid])
			fields[6] = "{0:.2f}".format(gfpkms[gid])
			fout.write("\t".join(fields[:7]) + "\n")

			gid += 1

	print("generate_expression_files is done.")



names, refs = load_ref(reference_fasta)
calc_background(names, refs, expression_file)
calc_foreground(names, refs, bam_file)
GCfactor = calcGCfactor(foreground, background)
generate_expression_files(refs, expression_file, gene_expression_file, out_expression_file, out_gene_expression_file)
