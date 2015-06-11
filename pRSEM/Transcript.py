__doc__="""

  peng 20131009

  Data structure copied from RSEM, made some name changes
"""


class Transcript:
  def __init__(self):
    self.transcript_id = None
    self.gene_id = None
    self.gene    = None
    self.transcript_group = None

    self.chrom   = None  ## RSEM Transcript's string seqname
    self.strand  = None
    self.length  = None
    self.exon_ranges = [];  ## RSEM Transcript's vector<Interval> structure
    self.gtf_attr = {}; ## RSEM Transcript's string left
    self.gtf_additional_info = None;

    self.start = None; ## genomic starting postion,
                       ## regardless of strand direction
                       ## always have a number smaller than self.end
    self.end = None;   ## genomic ending postion

   #self.expected_count = None; ## calculated by RSEM
   #self.tpm = None;            ## calculated by RSEM

    self.tss = None;   ## genomic coordinate of transcription starting site
    self.tes = None;   ## genomic coordinate of transcription ending site

   #self.nfrags = [];  ## number of stacked fragments at [TSS-500, end+500]
   #self.pos2nfrag = {}; ## nonzero number of stacked fragments for each
   #                     ## genomic position

    ## for Pol2 signal & peak
   #self.ave_signal_around_TSS  = None
   #self.ave_signal_around_body = None
   #self.ave_signal_around_TES  = None

   #self.has_peak_around_TSS  = None
   #self.has_peak_around_body = None
   #self.has_peak_around_TES  = None

   #self.peak_type = None ## defined based on if gene and grp has TSS peak

   #self.tss_region_contains_other_tss = None; ## True if [TSS-500, TSS+500]
                                               ## contains other transcript's
                                               ## TSS.
                                               ## False, otherwise
   #self.all_exons_overlap_with_other = None; ## True if all exons overlap with
   #                                          ## any other transcript
                                              ## False, otherwise

   #self.ave_mpp_over_exons = None ## average mappability over all exons

   #self.has_peak_around_TSS = None ## True, if has a peak overlap with
   #                                ## [TSS-500, TSS+500]; False, otherwise

    ## True if self's all exons+introns overlap another transcript's
    ## exons+introns
   #self.overlap_with_transcript_from_other_gene = None

    ## mappability
    self.ave_mpp_around_TSS = None  ## [TSS-flanking_width, TSS+flanking_width]
   #self.max_mpp_around_TSS = None

    self.ave_mpp_around_body = None ## (TSS+flanking_width, TES-flanking_width)
   #self.max_mpp_around_body = None

    self.ave_mpp_around_TES = None  ## [TES-flanking_width, TES+flanking_width]
   #self.max_mpp_around_TES = None


  def __str__(self):
    s = "%s\n%s\n%s\n%s %d\n" % (self.transcript_id, self.gene_id, self.chrom,
      self.strand, self.length);
    s += "%d" % len(self.exon_ranges);
    for (start, end) in self.exon_ranges:
      s += " %d %d" % (start, end);
    s += "\n";
    for key in self.gtf_attr.keys():
      for val in self.gtf_attr[key]:
        s += '%s "%s"; ' % (key, val);
    s = s.rstrip();
    return s;


  def constructFromRSEMTI(self, ti_lines):
    """
    construct Transcript from the 6 lines from RSEM .TI file
    """
    self.quicklyConstructFromRSEMTI(ti_lines);

    feature_words = ti_lines[5].rstrip(';').split(';');
    for feature_word in feature_words:
      feature_word.lstrip();
      (key, val) = feature_word.split();
      if not self.gtf_attr.has_key(key):
        self.gtf_attr[key] = [];
      self.gtf_attr[key].append(val.strip('"'));


  def quicklyConstructFromRSEMTI(self, ti_lines):
    """
    quickly construct Transcript from the 6 lines from RSEM .TI file, the last
    line won't be parsed.
    """
    self.transcript_id = ti_lines[0];
    self.gene_id = ti_lines[1];
    self.chrom = ti_lines[2];
    (self.strand, self.length) = ti_lines[3].split();
    self.length = int(self.length);
    words = ti_lines[4].split();
    for j in range(0, int(words[0])):
      start = int(words[j*2+1]);
      end = int(words[j*2+2]);
      self.exon_ranges.append( (start, end) );

    self.start = self.exon_ranges[0][0];
    self.end = self.exon_ranges[-1][-1];
    if self.strand == '+':
      self.tss = self.start
      self.tes = self.end
    elif self.strand == '-':
      self.tss = self.end
      self.tes = self.start
    self.gtf_additional_info = ti_lines[5];


  def defineTSSAndTES(self):
    """
    define TSS and TES
    """
    if (self.tss is None) or (self.tes is None):
      if self.strand == '+':
        self.tss = self.start;
        self.tes = self.end;
      elif self.strand == '-':
        self.tss = self.end;
        self.tes = self.start;


 #def ifTSSRegionContainsOtherTSS(self, other_transcripts, region_radius=None):
 #  """
 #  test if [TSS - region_radius, TSS + region_radius] contains other
 #  transcript's TSS
 #  self.tss_region_contains_other_tss = True if Yes
 #  self.tss_region_contains_other_tss = False if No
 #  """
 #  if region_radius is None:
 #    region_radius = 500

 #  if self.tss is None :
 #    self.defineTSSAndTES()

 #  self.tss_region_contains_other_tss = False
 #  for tr in other_transcripts:
 #    if (tr.chrom == self.chrom) and (tr.gene_id != self.gene_id):
 #      if tr.tss is None:
 #        tr.defineTSSAndTES()

 #      if abs(tr.tss - self.tss) < region_radius:
 #        self.tss_region_contains_other_tss = True
 #        break



 #def ifAllExonsOverlapWithOther(self, other_transcripts):
 #  """
 #  self.all_exons_overlap_with_other = True, if all exons in this transcripts
 #  are equal to or within other exons from the given list of transcripts,
 #  e.g. if half of self.exons overlap with transcript A's exons and the other
 #  half of self.exons overlap with transcript B's exons. Then it is True.

 #  self.all_exons_overlap_with_other = False, otherwise
 #  """
 #  import Util

 #  other_trs = filter(lambda tr: (tr.transcript_id != self.transcript_id) and
 #                     (tr.gene_id != self.gene_id) and
 #                     (tr.chrom == self.chrom),
 #                     other_transcripts)

 #  lol_ivs = [other_tr.exon_ranges for other_tr in other_trs]

 #  import itertools
 #  other_ivs = list(itertools.chain.from_iterable(lol_ivs))

 # #bools = map(lambda self_iv: Util.isSubsetOfIntervals(self_iv, other_ivs),
 # #            self.exon_ranges)
 # #self.all_exons_overlap_with_other = all(bools)

 #  self.all_exons_overlap_with_other = True
 #  for self_iv in self.exon_ranges:
 #    if not Util.isSubsetOfIntervals(self_iv, other_ivs):
 #      self.all_exons_overlap_with_other = False
 #      break


 #def ifOverlapWithTranscriptFromOtherGene(self, other_transcripts):
 #  """
 #  self.overlap_with_transcript_from_other_gene = True,
 #  if [self.start, self.end] is part or 100% overlap with another transcript
 #  False, otherwise
 #  """
 #  other_trs = filter(lambda tr: (tr.transcript_id != self.transcript_id) and
 #                     (tr.gene_id != self.gene_id) and
 #                     (tr.chrom == self.chrom),
 #                     other_transcripts)

 #  self.overlap_with_transcript_from_other_gene = \
 #    any( map(lambda tr: (self.start >= tr.start) and (self.end <= tr.end),
 #             other_trs))



 #def calculateExonMappability(self, genome=None):
 #  """
 #  calculate the average mappability for all exons in this transcript
 #  """
 #  import Util

 #  if genome is None:
 #    genome = 'hg19'

 #  sum_mpp = 0.0
 #  total_exon_len = 0
 #  for (start, end) in self.exon_ranges:
 #    if start == end: ## e.g. MCAM's 1st exon: chr11:119191799-119191799
 #      continue
 #    mean_mpp = Util.calculateMappability('mean', self.chrom, start, end,
 #                                         genome)
 #    sum_mpp += mean_mpp * ( end - start + 1)
 #    total_exon_len += end - start + 1

 #  self.ave_mpp_over_exons = sum_mpp*1.0/total_exon_len


  def calculateMappability(self, bin_bigwigsummary, fbigwig, width=500,
                           quiet=True):
    """
    calculate average mappability for a transcript's
    TSS region:  [TSS-width,   TSS+width],
    body region: [start+width+1, end-width-1],
    TES region:  [TES-width,   TES+width]

    if start+width+1 > end-width-1, then define body region as
      [end-width-1, start+width+1]

    assign the values for
    self.ave_mpp_around_TSS,  self.max_mpp_around_TSS
    self.ave_mpp_around_body, self.max_mpp_around_body
    self.ave_mpp_around_TES,  self.max_mpp_around_TES
    """
    import Util

    if (self.tss is None) or (self.tes is None):
      self.defineTSSAndTES()

    self.ave_mpp_around_TSS = Util.calculateMappability('mean', self.chrom,
                              self.tss - width, self.tss + width,
                              bin_bigwigsummary, fbigwig, quiet)

    if (self.start + width + 1) < (self.end - width - 1):
      self.ave_mpp_around_body = Util.calculateMappability('mean', self.chrom,
                                 self.start+width+1, self.end-width-1,
                                 bin_bigwigsummary, fbigwig, quiet)
    elif (self.start + width + 1) > (self.end - width - 1):
      self.ave_mpp_around_body = Util.calculateMappability('mean', self.chrom,
                                 self.end-width-1, self.start+width+1,
                                 bin_bigwigsummary, fbigwig, quiet)
    elif (self.start + width + 1) == (self.end - width - 1):
      self.ave_mpp_around_body = 1.0

    self.ave_mpp_around_TES = Util.calculateMappability('mean', self.chrom,
                              self.tes - width, self.tes + width,
                              bin_bigwigsummary, fbigwig, quiet)


 #def zeroChIPSeqSignal(self):
 #  self.ave_signal_around_TSS  = 0.0
 #  self.ave_signal_around_body = 0.0
 #  self.ave_signal_around_TES  = 0.0


 #def getChIPSeqSignal(self, stacked_fragment):
 #  """
 #  compute length-averaged ChIP-Seq signal for:
 #  [TSS-500, TSS+500], (TSS+500, TES-500), [TES-500, TES+500]

 #  signal is from StackedFragment defined in self.pos2nfrag
 #  """
 #  import numpy as np

 #  if len(self.gene.pos2nfrag) == 0:
 #    self.gene.assignStackedFragment(stacked_fragment)

 #  self.zeroChIPSeqSignal()

 #  if (self.end - self.start + 1 ) >= 1002:
 #    positions = []
 #    if self.strand == '+':
 #      positions = range(self.start-500, self.end+500+1)
 #    elif self.strand == '-':
 #      positions = range(self.end+500, self.start-500-1, -1)

 #    all_nfrags = np.zeros(len(positions))
 #    for (i, pos) in enumerate(positions):
 #      if pos in self.gene.pos2nfrag:
 #        all_nfrags[i] = self.gene.pos2nfrag[pos]

 #    self.ave_signal_around_TSS  = np.average(all_nfrags[:1001])
 #    if len(all_nfrags[1001:-1001]) > 0: ## in case of tr length == 1002
 #      self.ave_signal_around_body = np.average(all_nfrags[1001:-1001])
 #    self.ave_signal_around_TES  = np.average(all_nfrags[-1001:])



def readRSEMTI(fin):
  """
  read RSEM's .ti file, return a list of Transcripts objects
  """
  import Util

  lines = Util.readFile(fin);
  (ntranscripts, foo) = lines[0].split();
  ntranscripts = int(ntranscripts);
  transcripts = [];
  for i in range(0, ntranscripts):
    tr = Transcript();
    tr.constructFromRSEMTI(lines[i*6+1:i*6+7]);
    transcripts.append(tr);
   #if (i > 0) and (i % 20000 == 0):
   #  print "processed %d transcripts" % i;

  return transcripts;


def quicklyReadRSEMTI(fin):
  """
  read RSEM's .ti file without parsing the additional information line (the last
  line in a transcript's block

  return a list of Transcripts objects
  """
  import Util

  lines = Util.readFile(fin);
  (ntranscripts, foo) = lines[0].split();
  ntranscripts = int(ntranscripts);
  transcripts = [];
  for i in range(0, ntranscripts):
    tr = Transcript();
    tr.quicklyConstructFromRSEMTI(lines[i*6+1:i*6+7]);
    transcripts.append(tr);
    if (i > 0) and (i % 20000 == 0):
      print "processed %d transcripts" % i;

  return transcripts;


#def readPol2Features(fin):
# """
# read Pol2 feature files, create and return a list of transcripts
# """
# import Util

# trs = []
# lines = Util.readFile(fin)
# for line in lines[1:]:
#   words = line.split("\t")
#   tr = Transcript()
#   tr.transcript_id = words[0]
#   tr.gene_id       = words[1]
#   tr.chrom         = words[2]
#   tr.strand        = words[3]
#   tr.tss           = int(words[4])
#   tr.tes           = int(words[5])
#   tr.start         = int(words[6])
#   tr.end           = int(words[7])
#   tr.has_peak_around_TSS  = None
#   tr.has_peak_around_body = None
#   tr.has_peak_around_TES  = None

#   if words[8] == '1':
#     tr.has_peak_around_TSS = True
#   elif words[8] == '0':
#     tr.has_peak_around_TSS = False

#   if words[9] == '1':
#     tr.has_peak_around_body = True
#   elif words[9] == '0':
#     tr.has_peak_around_body = False

#   if words[10] == '1':
#     tr.has_peak_around_TES = True
#   elif words[10] == '0':
#     tr.has_peak_around_TES = False

#   tr.expected_count         = float(words[11])
#   tr.tpm                    = float(words[12])
#   tr.effective_length       = float(words[13])
#   tr.GC_content             = float(words[14])
#   tr.ave_signal_around_TSS  = float(words[15])
#   tr.ave_signal_around_body = float(words[16])
#   tr.ave_signal_around_TES  = float(words[17])

#   trs.append(tr)

# return trs

