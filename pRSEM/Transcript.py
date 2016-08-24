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

    self.tss = None;   ## genomic coordinate of transcription starting site
    self.tes = None;   ## genomic coordinate of transcription ending site

    ## mappability
    self.ave_mpp_around_TSS = None  ## [TSS-flanking_width, TSS+flanking_width]
    self.ave_mpp_around_body = None ## (TSS+flanking_width, TES-flanking_width)
    self.ave_mpp_around_TES = None  ## [TES-flanking_width, TES+flanking_width]


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
    self.transcript_id = ti_lines[0].split("\t")[0]
    self.gene_id = ti_lines[1].split("\t")[0]
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
    if (i > 0) and (i % 20000 == 0):
      print "processed %d transcripts" % i;

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

