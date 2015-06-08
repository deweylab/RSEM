__doc__="""

  pliu 20131024

  module for a group of Transcript objects.

  the initial purpose is for a group of transcripts that have TSS near each
  other
"""

class TranscriptGroup:
  def __init__(self):
    self.id = None
    self.chrom = None;
    self.gene_id = None;
    self.gene = None
    self.strand = None;
    self.ave_tss = None; ## averaged TSS position
    self.ave_tpm = None; ## averaged TPM
    self.transcripts = [];

    self.tss = None; ## strand: +, min(all_tss); strand: -, max(all_tss)
    self.tes = None; ## strand: +, max(all_tes); strand: -, min(all_tes)

   #self.nfrags = [];    ## pliu 20131126
                         ## need to implement it as a hash instead of array
    self.nfrag_start = None;
    self.nfrag_end = None;

    self.has_peak_around_TSS = None ## True, if has a peak overlap with any
                                    ## transcript's [TSS-500, TSS+500]
                                    ## False, otherwise

   #self.multi_tss_peaks = [] ## a list of peaks shared with other TSS groups
   #self.uni_tss_peaks = []   ## a list of peaks that do not share with others

 #def __str__(self):
 #  s = "Group: %s\t%s\t%s\t%.2f\t%.2f\n" % (self.chrom, self.strand,
 #    self.gene_id, self.ave_tss, self.ave_tpm);
 #  for tr in self.transcripts:
 #    s += "Transcript: %s\t%d\t%d\t%.2f\n" % (tr.transcript_id, tr.start,
 #      tr.end, tr.tpm);
 #    for pos in sorted(tr.pos2nfrag):
 #      if tr.pos2nfrag[pos] > 50:
 #        s += "%d\t%d\n" % (pos, tr.pos2nfrag[pos]);
 #  return s;


  def __str__(self):
    ss = [ "%s\t%s\t%s" % (self.gene_id, self.chrom, self.strand) ]
    if self.id is not None:
      ss += [ "\t%s" % self.id ]
    ss += [ "\n" ]
    for tr in self.transcripts:
      ss += [ "%s\t%d\n" % (tr.transcript_id, tr.tss) ]

    return ''.join(ss)


  def getTSSProfString(self):
    s = "%s\t%s\t%s\t%.2f\t" % (self.chrom, self.gene_id, self.strand,
      self.ave_tpm);
    for tr in self.transcripts:
      s += "%s," % tr.transcript_id;
    s = s.rstrip(',') + "\t";
    s += "%d\t%d\t" % (self.nfrag_start, self.nfrag_end);
   #for nfrag in self.nfrags:
   #  s += "%d\t" % nfrag;
    s = s.rstrip("\t") + "\n";
    return s;


  def calculateAverageTPM(self):
    self.ave_tpm = 0.0;
    for tr in self.transcripts:
      self.ave_tpm += tr.tpm;
    self.ave_tpm /= len(self.transcripts);


  def calculateAverageTSS(self):
    self.ave_tss = 0.0;
    for tr in self.transcripts:
      if tr.strand == '+':
        self.ave_tss += tr.start;
      elif tr.strand == '-':
        self.ave_tss += tr.end;
    self.ave_tss /= len(self.transcripts);


  def ifHasPeakAroundTSS(self, peaks):
    """
    peaks: a list of NarrowPeak objects
    True, if has a peak overlap with any transcript's [TSS-500, TSS+500]
    False, otherwise
    """
    import Util

    all_tss = [tr.tss for tr in self.transcripts]
    self.has_peak_around_TSS = Util.ifHasPeakOverlapWithInterval(peaks,
                                                               min(all_tss)-500,
                                                               max(all_tss)+500)

  def defineIfHasPeakAroundTSS(self):
    """
    a group is considered to have peak around TSS if any of its transcripts has
    peak around TSS
    """
    self.has_peak_around_TSS = False
    for tr in self.transcripts:
      if tr.has_peak_around_TSS:
        self.has_peak_around_TSS = True
        break


 #def calculateTSSAverageNFrags(self):
 #  """
 #  calculate average number of stacked fragments for [TSS-500, TSS+500];
 #  """
 #  import numpy as np;

 #  self.ave_tss = int(round(self.ave_tss));
 #  start = self.ave_tss - 500;
 #  end = self.ave_tss + 500;

 #  self.nfrags = np.zeros(1001, dtype=np.int32);
 #  for i in range(start, end+1):
 #    nfrag = 0.0;
 #    for tr in self.transcripts:
 #      if tr.pos2nfrag.has_key(i):
 #        nfrag += tr.pos2nfrag[i];
 #    nfrag = int(round(nfrag/len(self.transcripts)));
 #    self.nfrags[i-start] = nfrag;

 #  if self.strand == '+':
 #    self.nfrag_start = start;
 #    self.nfrag_end = end;
 #  elif self.strand == '-':
 #    self.nfrag_start = end;
 #    self.nfrag_end = start;
 #    self.nfrags = self.nfrags[::-1];


def constructByTSS(trs, cutoff=None):
  """
  given a list of transcripts from the same gene, classify them into TSS groups
  if a transcript's TSS is within a group's first member's
  [TSS-cutoff, TSS+cutoff], it will be classfied into that group
  """
  if cutoff is None:
    cutoff = 500

  tss_grp = TranscriptGroup()
  tss_grp.transcripts.append(trs[0])
  gene_tss_grps = [tss_grp]
  for tr in trs[1:]:
    is_assigned = False
    for grp in gene_tss_grps:
      if abs(tr.tss - grp.transcripts[0].tss) <= cutoff:
        grp.transcripts.append(tr)
        is_assigned = True
        break
    if not is_assigned:
      new_tss_grp = TranscriptGroup()
      new_tss_grp.transcripts.append(tr)
      gene_tss_grps.append(new_tss_grp)

  for (i, grp) in enumerate(gene_tss_grps):
    grp.chrom   = grp.transcripts[0].chrom
    grp.strand  = grp.transcripts[0].strand
    grp.gene_id = grp.transcripts[0].gene_id
    grp.id = grp.gene_id + '.' + str(i+1)
    grp.gene = grp.transcripts[0].gene
    grp.gene.transcript_groups.append(grp)
    for tr in grp.transcripts:
      tr.transcript_group = grp

  return gene_tss_grps



def constructByGeneAndTSS(transcripts, cutoff=None):
  """
  return a list of TranscriptGroup objects from input transcripts
  first classify transcripts by their gene_id, then build TSS groups for
  transcripts from the same gene
  """
  if cutoff is None:
    cutoff = 500

  gene_id_to_trs = {}
  for tr in transcripts:
    if tr.gene_id in gene_id_to_trs:
      gene_id_to_trs[tr.gene_id].append(tr)
    else:
      gene_id_to_trs[tr.gene_id] = [tr]

  tss_grps = []
  for (gene_id, trs) in gene_id_to_trs.items():
    tss_grps += constructByTSS(trs, cutoff)

  return tss_grps
