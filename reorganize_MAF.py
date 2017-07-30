#!/usr/bin/env python

from sys import argv, exit

if len(argv) != 3:
	print("Usage: python reorganize.py input.mutation.txt output.maf")
	exit(-1)

class MAFLine:
	def __init__(self, pos, strlen, fields):

		self.pos = pos
		self.strlen = strlen
		self.fields = fields


def outputChr(fout, ins_list, del_list):

	m = len(ins_list)
	n = len(del_list)

	if m + n == 0:
		return

	i = 0
	for del_line in del_list:
		while i < m and ins_list[i].pos < del_line.pos - 1:
			fout.write("\t".join(ins_list[i].fields) + "\n")
			i += 1
		if i < m and ins_list[i].pos < del_line.pos + del_line.strlen:
			# print("ins_pos = {}, del_pos = {}".format(ins_list[i].pos, del_line.pos))
			
			ins_line = ins_list[i]
			ins_line.pos = del_line.pos - 1
			ins_line.fields[5] = str(ins_line.pos)
			ins_line.fields[6] = str(ins_line.pos + 1)
			i += 1
			while i < m and ins_list[i].pos < del_line.pos + del_line.strlen:
				ins_line.strlen += ins_list[i].strlen
				ins_line.fields[11] += ins_list[i].fields[11]
				ins_line.fields[12] += ins_list[i].fields[12]
				i += 1
			fout.write("\t".join(ins_line.fields) + "\n")
		fout.write("\t".join(del_line.fields) + "\n")

	while i < m:
		fout.write("\t".join(ins_list[i].fields) + "\n")
		i += 1
	

with open(argv[1]) as fin, open(argv[2], "w") as fout:
	header = next(fin)
	fout.write(header)

	cur_chrn = ""
	seen_chrs = set()
	ins_list = []
	del_list = []

	for line in fin:
		fields = line.strip().split('\t')

		chrn = fields[4]
		pos = int(fields[5])
		event = fields[9]
		strlen = len(fields[11]) if event == "INS" else len(fields[10])

		if chrn != cur_chrn:
			assert chrn not in seen_chrs

			outputChr(fout, ins_list, del_list)
			cur_chrn = chrn
			seen_chrs.add(chrn)
			ins_list = []
			del_list = []

		if event == "INS":
			assert len(fields[11]) == len(fields[12])

			if len(ins_list) > 0 and ins_list[-1].pos == pos:
				ins_list[-1].strlen += strlen
				ins_list[-1].fields[11] += fields[11]
				ins_list[-1].fields[12] += fields[12]
			else:
				ins_list.append(MAFLine(pos, strlen, fields))
		else:
			assert event == "DEL"
			if len(del_list) > 0 and del_list[-1].pos + del_list[-1].strlen >= pos:
				newlen = pos + strlen - del_list[-1].pos
				if newlen > del_list[-1].strlen:
					del_list[-1].fields[10] += fields[10][del_list[-1].pos + del_list[-1].strlen - pos : ]
					del_list[-1].strlen = newlen
			else:
				del_list.append(MAFLine(pos, strlen, fields))

	outputChr(fout, ins_list, del_list)
