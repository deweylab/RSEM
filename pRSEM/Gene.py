__doc__="""

  peng 20131002

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
   #self.narrow_peaks = [];
   #self.pos2nfrag = {};
    self.transcript_tss_groups = []; ## a list of list of transcripts having TSS
                                     ## within user-specified distance
    self.transcript_groups = [] ## a list of TranscriptionGroup objects

   #self.has_peak_around_TSS = None; ## True if has peak overlaped with
   #                                 ## [TSS-500, TSS+500]; False otherwise
   #                                 ## it is the same as
   #                                 ## self.has_pol2_peak_around_TSS
   #                                 ## but in a shorter name
   #self.has_peak_around_body = None; ## True if has peak overlaped with
   #                                  ## [TSS+500, TES-500]; False otherwise
   #self.has_peak_around_TES  = None; ## True if has peak overlaped with
   #                                  ## [TES-500, TSS+500]; False otherwise

   #self.peak_type = None ## 'all':   all grps have peaks
   #                      ## 'none':  no grp has peak
   #                      ## 'mixed': some grps has peaks, the others do not

   #self.tpm = None;
   #self.expected_count = None;   ## from RSEM's .gene.results
   #self.effective_length = None; ## from RSEM's .gene.results

   #self.GC_content = None; ## percentage of G or C in a sequence

   #self.cluster_id = None; ## from kmean cluster

   #self.pol2_bins = []; ## 300 bins for pol2 profile,
   #                     ## 1st 100 for [TSS - 500, TSS + 500]
   #                     ## 2nd 100 for (TSS + 500, TES - 500)
   #                     ## 3rd 100 for [TES - 500, TES + 500]

   #self.ave_signal_around_TSS  = None ## len-averaged signal [TSS-500, TSS+500]
   #self.ave_signal_around_body = None ## len-averaged signal [TSS+500, TES-500]
   #self.ave_signal_around_TES  = None ## len-averaged signal [TES-500, TES+500]

    ## mappability
   #self.ave_mpp_around_TSS = None;  ## [TSS-500, TSS+500]
   #self.max_mpp_around_TSS = None;

   #self.ave_mpp_around_body = None; ## (TSS+500, TES-500)
   #self.max_mpp_around_body = None;

   #self.ave_mpp_around_TES = None;  ## [TES-500, TES+500]
   #self.max_mpp_around_TES = None;

   #self.ave_mpp_over_exons = None ## average mpp over all exons
   #self.ave_mpp_over_transcripts = None ## average mpp over all transcripts

   #self.exon_ranges = [] ## in RSEM Transcripts's vector<Interval> format

   #self.has_mixed_peak = None ## True, some of its TSS groups have peak and
   #                           ## some do not
   #                           ## False, otherwise

   #self.has_TSS_close_to_other_gene = None ## True, if any TSS within 500 bp
   #                                        ## to TSS of other gene
   #                                        ## False, otherwise


  def __str__(self):
    s = "%s" % self.gene_id;
    return s;


 #def assignStackedFragment(self, stacked_fragment):
 #  """
 #  from input list of StackedFragment objects, find out numbers of stacked
 #  fragment for [TSS-500, TES+500]
 #  """
 #  self.pos2nfrag = {}
 #  wind_start = self.start - 500
 #  wind_end = self.end + 500
 #  for pos in range(wind_start, wind_end+1):
 #    if pos in stacked_fragment.pos2nfrag:
 #      self.pos2nfrag[pos] = stacked_fragment.pos2nfrag[pos]


 #def initChIPSeqSignal(self):
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

 #  if len(self.pos2nfrag) == 0:
 #    self.assignStackedFragment(stacked_fragment)

 #  self.initChIPSeqSignal()

 #  if (self.end - self.start + 1 ) >= 1002:
 #    positions = []
 #    if self.strand == '+':
 #      positions = range(self.start-500, self.end+500+1)
 #    elif self.strand == '-':
 #      positions = range(self.end+500, self.start-500-1, -1)

 #    all_nfrags = np.zeros(len(positions))
 #    for (i, pos) in enumerate(positions):
 #      if pos in self.pos2nfrag:
 #        all_nfrags[i] = self.pos2nfrag[pos]

 #    self.ave_signal_around_TSS  = np.average(all_nfrags[:1001])
 #    self.ave_signal_around_body = np.average(all_nfrags[1001:-1001])
 #    self.ave_signal_around_TES  = np.average(all_nfrags[-1001:])



 #def normalizeStackedFragment(self):
 #  """
 #  get the average of number of stacked fragments for each bin,
 #  save it in self.pol2_bins in the direction of TSS to TES

 #  [TSS-500, TSS+500] 100 bins
 #  (TSS+500, TES-500) 100 bins
 #  [TES-500, TES+500] 100 bins
 #  """
 #  import numpy as np;

 #  self.pol2_bins = np.zeros(300);

 #  wind_start = self.start - 500;
 #  wind_end = self.end + 500;

 #  length = self.end - self.start + 1;
 #  if length >= 1102:
 #    positions = [];
 #    if self.strand == '+':
 #      positions = range(wind_start, wind_end+1);
 #    elif self.strand == '-':
 #      positions = range(wind_end, wind_start-1, -1);

 #    all_nfrags = np.zeros(len(positions));
 #    for (i, pos) in enumerate(positions):
 #      if pos in self.pos2nfrag:
 #        all_nfrags[i] = self.pos2nfrag[pos];

 #    loc_types = ['tss', 'mid', 'end'];
 #    nfrags = {
 #      'tss': all_nfrags[:1001],
 #      'mid': all_nfrags[1001:-1001],
 #      'end': all_nfrags[-1001:],
 #    }

 #    for (i, loc_type) in enumerate(loc_types):
 #      bins = np.array_split(nfrags[loc_type], 100);
 #      for (j, bin) in enumerate(bins):
 #       #self.pol2_bins.append(np.average(bin)); ## average by the # of bps in
 #       #                                        ## each bin
 #        self.pol2_bins[i*100+j] = np.average(bin); ## averaged by the # of bps
 #                                                   ## in each bin
 #        if np.isnan(np.average(bin)):
 #          print self.end - self.start + 1, len(all_nfrags), all_nfrags;
 #          print loc_type, nfrags[loc_type];
 #          print bin;
 #          print np.average(bin);



 #def calculateMappability(self):
 #  """
 #  calculate the average and max mappability for
 #  [TSS-500, TSS+500], (TSS+500, TES-500), [TES-500, TES+500]

 #  assign the values for
 #  self.ave_mpp_around_TSS,  self.max_mpp_around_TSS
 #  self.ave_mpp_around_body, self.max_mpp_around_body
 #  self.ave_mpp_around_TES,  self.max_mpp_around_TES
 #  """
 #  import Util;

 #  if (self.tss is None) or (self.tes is None):
 #    self.getStartEndTSSTESFromTranscripts();

 #  self.ave_mpp_around_TSS  = Util.calculateMappability('mean', self.chrom,
 #    str(self.tss - 500), str(self.tss + 500) );
 # #self.max_mpp_around_TSS  = Util.calculateMappability('max',  self.chrom,
 # #  str(self.tss - 500), str(self.tss + 500) );

 #  self.ave_mpp_around_body = Util.calculateMappability('mean', self.chrom,
 #    str(self.start + 500), str(self.end - 500) );
 # #self.max_mpp_around_body = Util.calculateMappability('max',  self.chrom,
 # #  str(self.start + 500), str(self.end - 500) );

 #  self.ave_mpp_around_TES  = Util.calculateMappability('mean', self.chrom,
 #    str(self.tes - 500), str(self.tes + 500) );
 # #self.max_mpp_around_TES  = Util.calculateMappability('max',  self.chrom,
 # #  str(self.tes - 500), str(self.tes + 500) );


 #def calculateExonMappability(self):
 #  """
 #  calculate the average mappability for all exons in this gene
 #  """
 #  import Util

 #  if len(self.exon_ranges) == 0:
 #    self.getExonRanges()

 #  sum_mpp = 0.0
 #  total_exon_len = 0
 #  for (start, end) in self.exon_ranges:
 #    if start == end: ## e.g. MCAM's 1st exon: chr11:119191799-119191799
 #      continue
 #    mean_mpp = Util.calculateMappability('mean', self.chrom, start, end)
 #    sum_mpp += mean_mpp * ( end - start + 1)
 #    total_exon_len += end - start + 1

 #  self.ave_mpp_over_exons = sum_mpp*1.0/total_exon_len


 #def calculateAverageMappabilityOverTranscripts(self):
 #  """
 #  calculate the mean mappability over all transcripts in this gene
 #  """
 #  import numpy as np

 #  map(lambda tr: tr.calculateExonMappability(), self.transcripts)
 #  mpps = [tr.ave_mpp_over_exons for tr in self.transcripts]
 #  self.ave_mpp_over_transcripts = np.mean(mpps)



  ## should be moved to TranscriptGroup file
  def groupTranscriptsByTSS(self):
    """
    put transcripts that have TSS within certain distance into a group
    """
   #import Param;

   #cutoff = 100; ## TSS within 100 bp
   #cutoff = Param.TSS_GROUP_CUTOFF; ## cutoff for grouping TSS
    cutoff = 500 ## cutoff for grouping TSS

    group = [self.transcripts[0]]
    self.transcript_tss_groups.append(group)
    for tr in self.transcripts[1:]:
      is_assigned = False
      for grp in self.transcript_tss_groups:
        if (self.strand == '+') and (abs(tr.start - grp[0].start)<=cutoff):
          grp.append(tr)
          is_assigned = True;
        elif (self.strand == '-') and (abs(tr.end - grp[0].end)<=cutoff):
          grp.append(tr)
          is_assigned = True

        if is_assigned:
          break

      if not is_assigned:
        self.transcript_tss_groups.append([tr])


  ## should be moved to TranscriptGroup file
  def constructTranscriptGroups(self):
    """
    construct a list of TranscriptGroup objects
    """
    import TranscriptGroup

    if len(self.transcript_tss_groups) == 0:
      self.groupTranscriptsByTSS()

    for transcripts in self.transcript_tss_groups:
      grp = TranscriptGroup.TranscriptGroup()
      grp.chrom = self.chrom
      grp.gene_id = self.gene_id
      grp.strand = self.strand
      grp.transcripts = transcripts
      self.transcript_groups.append(grp)



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



 #def ifAllTranscriptsHaveTheSameStartEnd(self):
 #  """
 #  return True if all transcripts shared the same starting and ending positions
 #  return False otherwise
 #  """
 #  starts = [tr.start for tr in self.transcripts];
 #  ends = [tr.end for tr in self.transcripts];
 #  all_the_same = bool( (starts.count(starts[0]) == len(starts)) and \
 #                       (ends.count(ends[0])== len(ends)) );

 #  return all_the_same;


 #def getExonRanges(self):
 #  """
 #  if one position is an exon in any transcript, it will be an exon for this
 #  gene
 #  """
 #  dict_exon_positions = {}
 #  for tr in self.transcripts:
 #    for (start, end) in tr.exon_ranges:
 #      for i in range(start, end+1):
 #        dict_exon_positions.setdefault(i, True)

 #  exon_positions = sorted(dict_exon_positions.keys())

 #  position_arrays = [ [exon_positions[0]] ]
 #  for i in range(1, len(exon_positions)):
 #    pos = exon_positions[i]
 #    if (pos - position_arrays[-1][-1]) == 1:
 #      position_arrays[-1].append(pos)
 #    else:
 #      position_arrays.append([pos])

 #  for position_array in position_arrays:
 #    self.exon_ranges.append((position_array[0], position_array[-1]))


 #def ifHasTSSCloseToOtherGene(self, genes, dist_cutoff=None):
 #  """
 #  define self.has_TSS_close_to_other_gene = True, this any TSS from this gene
 #  is within dist_cutoff (default = 500bp) to TSS from other gene
 #  False, otherwise
 #  """
 #  if dist_cutoff is None:
 #    dist_cutoff = 500

 #  self.has_TSS_close_to_other_gene = False

 #  for gene in genes:
 #    if (gene.chrom == self.chrom) and (gene.gene_id != self.gene_id) and \
 #      (not self.has_TSS_close_to_other_gene):
 #      for other_tr in gene.transcripts:
 #        if self.has_TSS_close_to_other_gene:
 #          break
 #        for tr in self.transcripts:
 #          if abs(other_tr.tss - tr.tss) <= dist_cutoff:
 #           #print tr.transcript_id, other_tr.transcript_id
 #            self.has_TSS_close_to_other_gene = True
 #            break

 # #print self.gene_id, "\t", self.has_TSS_close_to_other_gene


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



#def constructGenesFromRSEMGeneResults(rsem_gene_results):
# """
# return a list of Gene objects constructed from a list of RSEMGeneResult
# """
# genes = [];
# for result in rsem_gene_results:
#   gene = Gene();
#   gene.gene_id = result.gene_id;
#   gene.tpm = result.tpm;
#   gene.expected_count = result.expected_count;
#   gene.effective_length = result.effective_length;
#   gene.rsem_result = result;
#   genes.append(gene);

# return genes;


def buildTrainingSet(genes, prm):
  """
  write training set in file Param.ftrainingg_tr
  """
  import Util

  ogot_genes = filter(lambda g: len(g.transcripts) == 1 and
                                (g.end - g.start + 1) >=
                                prm.TRAINING_GENE_MIN_LEN, genes)

  trs = [tr for g in ogot_genes for tr in g.transcripts]

  trid2mpps = Util.runMPOverAList(prm.num_threads, calTSSBodyTESMappability,
                                  [trs, prm])

  with open(prm.falltrcrd, 'w') as f_fout:
    f_fout.write("geneid\ttrid\tchrom\tstrand\tstart\tend\t")
    f_fout.write("tss_mpp\tbody_mpp\ttes_mpp\n")
    for gene in genes:
      for tr in gene.transcripts:
        f_fout.write("%s\t%s\t%s\t%s\t%d\t%d\t" % ( tr.gene_id,
                     tr.transcript_id, tr.chrom, tr.strand, tr.start, tr.end))
        if tr.transcript_id in trid2mpps:
          mpps = trid2mpps[tr.transcript_id]
          f_fout.write("%4.2f\t%4.2f\t%4.2f\n" % mpps)
        else:
          f_fout.write("NA\tNA\tNA\n")

  Util.runCommand('/bin/env', 'Rscript', prm.rnaseq_rscript, 'selTrainingTr',
                  prm.prsem_rlib_dir, prm.falltrcrd,
                  prm.TRAINING_MIN_MAPPABILITY, prm.FLANKING_WIDTH,
                  prm.ftraining_tr, quiet=prm.quiet)


def calTSSBodyTESMappability(trs, prm, out_q):
  """
  calculate average mappability around TSS, body, and TES for all transcripts of
  given list of genes

  save results in transcript's attribute
  """
  outdict = {}
  for tr in trs:
    tr.calculateMappability(prm.bigwigsummary_bin, prm.mappability_bigwig_file,
                            prm.FLANKING_WIDTH)
    outdict[tr.transcript_id] = (tr.ave_mpp_around_TSS, tr.ave_mpp_around_body,
                                 tr.ave_mpp_around_TES)
  out_q.put(outdict)
