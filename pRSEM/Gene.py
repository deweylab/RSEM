__doc__="""

  pliu 20131002

  module for gene
"""


class Gene:
  def __init__(self):
    self.gene_id = None;
   #self.rsem_result = None;

    self.chrom = None;
    self.strand = None;
    self.tss = None;   ## transcription starting sites, strand ori considered
    self.tes = None;   ## transcription ending sites, strand ori considered
    self.start = None; ## genomic starting position regardless of strand
                       ## direction, always have a number smaller than self.end
    self.end = None;   ## genomic ending position


    self.transcripts = [];
    self.gtfs = [];
    self.transcript_tss_groups = []; ## a list of list of transcripts having TSS
                                     ## within user-specified distance
    self.transcript_groups = [] ## a list of TranscriptionGroup objects


  def __str__(self):
    s = "%s" % self.gene_id;
    return s;


  ## should be moved to TranscriptGroup file
 #def groupTranscriptsByTSS(self):
 #  """
 #  put transcripts that have TSS within certain distance into a group
 #  """
 # #import Param;

 # #cutoff = 100; ## TSS within 100 bp
 # #cutoff = Param.TSS_GROUP_CUTOFF; ## cutoff for grouping TSS
 #  cutoff = 500 ## cutoff for grouping TSS

 #  group = [self.transcripts[0]]
 #  self.transcript_tss_groups.append(group)
 #  for tr in self.transcripts[1:]:
 #    is_assigned = False
 #    for grp in self.transcript_tss_groups:
 #      if (self.strand == '+') and (abs(tr.start - grp[0].start)<=cutoff):
 #        grp.append(tr)
 #        is_assigned = True;
 #      elif (self.strand == '-') and (abs(tr.end - grp[0].end)<=cutoff):
 #        grp.append(tr)
 #        is_assigned = True

 #      if is_assigned:
 #        break

 #    if not is_assigned:
 #      self.transcript_tss_groups.append([tr])


  ## should be moved to TranscriptGroup file
 #def constructTranscriptGroups(self):
 #  """
 #  construct a list of TranscriptGroup objects
 #  """
 #  import TranscriptGroup

 #  if len(self.transcript_tss_groups) == 0:
 #    self.groupTranscriptsByTSS()

 #  for transcripts in self.transcript_tss_groups:
 #    grp = TranscriptGroup.TranscriptGroup()
 #    grp.chrom = self.chrom
 #    grp.gene_id = self.gene_id
 #    grp.strand = self.strand
 #    grp.transcripts = transcripts
 #    self.transcript_groups.append(grp)



  def getStartEndTSSTESFromTranscripts(self):
    """
    define start and end from gene's transcripts
    start = min{all starts for transcripts};
    end = max{all ends for transcripts};
    """
    starts = [tr.start for tr in self.transcripts];
    ends = [tr.end for tr in self.transcripts];
    self.start = min(starts);
    self.end = max(ends);

    if self.strand == '+':
      self.tss = self.start;
      self.tes = self.end;
    elif self.strand == '-':
      self.tss = self.end;
      self.tes = self.start;


  def definePeakTypeByTranscriptGroups(self):
    """
    all:  all its transcript groups have peaks
    none: none of its transcript groups has peak
    mixed: some of its transcript groups have peaks, the others do not
    """
    has_tss_peaks = [grp.has_peak_around_TSS for grp in self.transcript_groups]
    if all(has_tss_peaks): ## all groups have peaks
      self.peak_type = 'all'
    else:
      if any(has_tss_peaks): ## some groups have peaks, the others not
        self.peak_type = 'mixed'
      else:                  ## no group has peak
        self.peak_type = 'no'


def constructGenesFromTranscripts(transcripts):
  """
  return a list of genes constructed from input transcripts
  """
  genes = []
  gene_dict_id = {}
  for tr in transcripts:
    if gene_dict_id.has_key(tr.gene_id):
      gene_dict_id[tr.gene_id].transcripts.append(tr)
      tr.gene = gene_dict_id[tr.gene_id]
    else:
      gene = Gene()
      gene.gene_id = tr.gene_id
      gene.chrom = tr.chrom
      gene.strand = tr.strand
      gene.transcripts.append(tr)
      genes.append(gene)
      gene_dict_id[tr.gene_id] = gene
      tr.gene = gene

  map(lambda gene: gene.getStartEndTSSTESFromTranscripts(), genes);

  return genes;

