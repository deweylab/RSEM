__doc__="""

  pliu 20150511

  python module for a ChIP-seq experiment that contains
  replicates of ChIP-seq data for target and/or control
"""

import ChIPSeqReplicate
import File
import Util


class ChIPSeqExperiment:
  def __init__(self):
    self.param           = None ## reference to input parameters
    self.reps            = []   ## list of ChIPSeqReplciate object
    self.is_control      = None ## if is control
    self.pooled_tagalign = None ## File object of pooled tagAlign


  @classmethod
  def initFromParam(cls, param, is_control, param_attr):
    cse = cls()
    cse.param = param
    cse.is_control = is_control
    ftgts = getattr(param, param_attr).split(',')
    cse.reps = [ ChIPSeqReplicate.initFromFastqFile(ffq) for ffq in ftgts ]
    for (i, rep) in enumerate(cse.reps):
      rep.param = param
      rep.index = i+1
      rep.chipseqexp = cse
      tgt_fta = "%s/%s.tagAlign.gz" % (param.temp_dir, rep.name)
      rep.tagalign = File.initFromFullFileName(tgt_fta)
    return cse


  def getFastqEncoding(self):
    nthr = self.param.num_threads
    fin = ','.join([ f.fastq.fullname for f in self.reps])
    if self.is_control:
      fenc = self.param.imd_name + '_prsem.chipseq_control_encoding'
    else:
      fenc = self.param.imd_name + '_prsem.chipseq_target_encoding'

    Util.runCommand(self.param.chipseq_rscript, 'guessFqEncoding',
                    nthr, fin, fenc, self.param.prsem_rlib_dir)

    with open(fenc, 'r') as f_fenc:
      next(f_fenc)
      file2enc = dict([ line.rstrip("\n").split("\t") for line in f_fenc ])

    for f in self.reps:
      f.encoding = file2enc[f.fastq.fullname]


  def alignReadByBowtie(self):
    import os
    if self.param.num_threads > 4:
      nthr_bowtie = self.param.num_threads - 4
    else:
      nthr_bowtie = 1

    bowtie_ref_name = "%s_prsem" % self.param.ref_name
    for rep in self.reps:
      cmd_cat = Util.getCatCommand(rep.fastq.is_gz)

      ## many pipes, have to use os.system
      cmds = [ "%s %s |" % (cmd_cat, rep.fastq.fullname) ] + \
             [ "%s -q -v 2 -a --best --strata -m 1 %s -S -p %d %s - |" % (
               self.param.bowtie_bin_for_chipseq, rep.encoding, nthr_bowtie,
               bowtie_ref_name ) ] + \
             [ "%s view -S -b -F 1548 - |" % self.param.samtools_bin ] + \
             [ "%s bamtobed -i stdin |" % (
               self.param.bedtools_bin_for_chipseq ) ] + \
             [ """awk 'BEGIN{FS="\\t";OFS="\\t"}{$4="N"; print $0}' |""" ] + \
             [ "gzip -c > %s " % rep.tagalign.fullname ]

      cmd = ' '.join(cmds)
      print cmd, "\n";
      os.system(cmd)


  def poolTagAlign(self):
    import os
    if self.is_control:
      frep0 = self.param.temp_dir + 'controlRep0.tagAlign.gz'
    else:
      frep0 = self.param.temp_dir + 'targetRep0.tagAlign.gz'

    if os.path.exists(frep0):
      os.remove(frep0)

    self.pooled_tagalign = File.initFromFullFileName(frep0)
    for rep in self.reps:
      cat_cmd = Util.getCatCommand(rep.fastq.is_gz)
      cmd = "%s %s | gzip -c >> %s" % (cat_cmd, rep.tagalign.fullname, frep0)
      print cmd, "\n";
      os.system(cmd)


  def runSPP(self, fctrl_tagalign):
    """
    in principle, this function is only for ChIP-seq target experiment
    should make target and control inherit from ChIPSeqExperiment, do it later
    """
    import sys
    if self.is_control:
      sys.exit( "ChIPSeqExperiment::runSPP() cann't be applied to control" )

    tgt_tagaligns = [self.pooled_tagalign] + [rep.tagalign for rep in self.reps]
    prm = self.param

    ## need to check and install spp ##
    Util.runCommand(prm.chipseq_rscript, 'checkInstallSpp', prm.spp_tgz,
                    prm.prsem_rlib_dir)
    exit(-1)

    for tgt_tagalign in tgt_tagaligns:
      fout = "%sphantom_%s.tab" % (prm.temp_dir, tgt_tagalign.basename)
      Util.runCommand('/bin/env', 'Rscript', prm.spp_script,
                      "-c=%s"      % tgt_tagalign.fullname,
                      "-i=%s"      % fctrl_tagalign,
                      "-npeak=%d"  % prm.N_PEAK,
                      prm.PEAK_TYPE,
                      '-savp',
                      "-x=%s"      % prm.EXCLUSION_ZONE,
                      '-rf',
                      "-odir=%s"   % prm.temp_dir,
                      "-p=%d"      % prm.num_threads,
                      "-tmpdir=%s" % prm.temp_dir,
                      "-out=%s"    % fout )



def initFromParam(param, typ):
  if typ.lower() == 'target':
    is_ctrl = False
    param_attr = 'chipseq_target_read_files'
  elif typ.lower() in [ 'control', 'input' ]:
    is_ctrl = True
    param_attr = 'chipseq_control_read_files'

  return ChIPSeqExperiment.initFromParam(param, is_ctrl, param_attr)
