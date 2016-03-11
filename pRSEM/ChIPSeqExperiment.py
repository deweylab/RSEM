__doc__="""

  pliu 20150511

  python module for a ChIP-seq experiment that contains
  replicates of ChIP-seq data for target and/or control
"""

import os
import sys
import multiprocessing as mp

import Util


class ChIPSeqExperiment:
  def __init__(self):
    self.param           = None ## reference to input parameters
    self.reps            = []   ## list of ChIPSeqReplcate object
    self.is_control      = None ## if is control
    self.pooled_tagalign = None ## File obj of pooled tagAlign
    self.peaks           = None ## File obj of targetRep0_VS_controlRep0 peaks
    self.final_peaks     = None ## File obj of final peaks


  @classmethod
  def initFromParam(cls, param, is_control, param_attr):
    import ChIPSeqReplicate
    import File
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

    if cse.is_control:
      frep0 = param.fchipseq_control_signals
    else:
      frep0 = param.fchipseq_target_signals
    cse.pooled_tagalign = File.initFromFullFileName(frep0)

    cse.peaks = File.initFromFullFileName(param.fall_chipseq_peaks)

    cse.final_peaks = File.initFromFullFileName(param.fidr_chipseq_peaks)

    return cse


  def getFastqEncoding(self):
    nthr = self.param.num_threads
    fin = ','.join([ f.fastq.fullname for f in self.reps])
    if self.is_control:
      fenc = self.param.imd_name + '_prsem.chipseq_control_encoding'
    else:
      fenc = self.param.imd_name + '_prsem.chipseq_target_encoding'

    Util.runCommand('/bin/env', 'Rscript', self.param.chipseq_rscript,
                    'guessFqEncoding', nthr, fin, fenc,
                    self.param.prsem_rlib_dir, quiet=self.param.quiet )

    if not os.path.exists(fenc):
      sys.exit("Failed to generate file: %s\n" % fenc)

    with open(fenc, 'r') as f_fenc:
      next(f_fenc)
      file2enc = dict([ line.rstrip("\n").split("\t") for line in f_fenc ])

    for f in self.reps:
      f.encoding = file2enc[f.fastq.fullname]


  def alignReadByBowtie(self):
    ## running zat, filterSam2Bed and gzip takes about 1 thread
    if self.param.num_threads > 2:
      nthr_bowtie = self.param.num_threads - 1
    else:
      nthr_bowtie = 1

    bowtie_ref_name = "%s_prsem" % self.param.ref_name
    for rep in self.reps:
      cmd_cat = Util.getCatCommand(rep.fastq.is_gz)

      if not os.path.exists(rep.fastq.fullname):
        sys.exit("File not found: %s\n" % rep.fastq.fullname)

      s_quiet = None
      if self.param.quiet:
        s_quiet = ' --quiet '
      else:
        s_quiet = ''

      ## many pipes, have to use os.system
      cmds = [ "%s %s |" % (cmd_cat, rep.fastq.fullname) ] + \
             [ "%s/bowtie " % self.param.bowtie_path ]  + \
             [ "%s -q -v 2 -a --best --strata -m 1 %s -S -p %d %s - | " % (
               s_quiet, rep.encoding, nthr_bowtie, bowtie_ref_name ) ] + \
             [ "%s - | " % self.param.filterSam2Bed ] + \
             [ "gzip -c > %s " % rep.tagalign.fullname ]

      cmd = ' '.join(cmds)

      ## use all threads to align ChIP-seq reads sequentially
      Util.runCommand(cmd, quiet=self.param.quiet)

      if not os.path.exists(rep.tagalign.fullname):
        sys.exit("failed to generate file: %s\n" % rep.tagalign.fullname)



  def poolTagAlign(self):
    frep0 = self.pooled_tagalign.fullname
    if os.path.exists(frep0):
      os.remove(frep0)
    for rep in self.reps:
      cat_cmd = Util.getCatCommand(rep.fastq.is_gz)
      if not os.path.exists(rep.tagalign.fullname):
        sys.exit("File not found: %s\n" % rep.tagalign.fullname)

      cmd = "%s %s | gzip -c >> %s" % (cat_cmd, rep.tagalign.fullname, frep0)
      Util.runCommand(cmd, quiet=self.param.quiet)

      if not os.path.exists(frep0):
        sys.exit("Failed to generate file: %s\n" % frep0)


  def callPeaksBySPP(self, ctrl_tagalign):
    """
    in principle, this function is only for ChIP-seq target experiment
    should make target and control inherit from ChIPSeqExperiment, will do
    """
    if self.is_control:
      sys.exit( "ChIPSeqExperiment::runSPP() cann't be applied to control" )

    tgt_tagaligns = [self.pooled_tagalign] + [rep.tagalign for rep in self.reps]
    prm = self.param

    ## need to add pRSEM's R_LIBS path so that run_spp.R can load spp library
    if 'R_LIBS' in os.environ:
      os.environ['R_LIBS'] = "%s:%s" % (os.environ['R_LIBS'],
                                        prm.prsem_rlib_dir)
    else:
      os.environ['R_LIBS'] = prm.prsem_rlib_dir

    nthr = prm.num_threads/len(tgt_tagaligns)
    fctrl_tagalign = ctrl_tagalign.fullname
    procs = [ mp.Process(target=runSPP, args=(tgt_tagalign, fctrl_tagalign,
                         prm, nthr)) for tgt_tagalign in tgt_tagaligns ]
    for p in procs:
      p.start()
    for p in procs:
      p.join()


  def getPeaksByIDR(self, ctrl_tagalign):
    """
    in principle, this function is only for ChIP-seq target experiment
    should make target and control inherit from ChIPSeqExperiment, will do
    """
    import gzip
    import itertools
    if self.is_control:
      sys.exit( "ChIPSeqExperiment::runSPP() can't be applied to control" )

    procs = []
    out_q = mp.Queue()
    prm = self.param
    for (repa, repb) in itertools.combinations(self.reps, 2):
      fpeaka = prm.temp_dir + repa.tagalign.filename_sans_ext + '_VS_' + \
               ctrl_tagalign.filename_sans_ext + '.regionPeak.gz'
      fpeakb = prm.temp_dir + repb.tagalign.filename_sans_ext + '_VS_' + \
               ctrl_tagalign.filename_sans_ext + '.regionPeak.gz'
      if not os.path.exists(fpeaka):
        sys.exit("File not found: %s\n" % fpeaka)
      if not os.path.exists(fpeakb):
        sys.exit("File not found: %s\n" % fpeakb)

      idr_prefix = prm.temp_dir + 'idr_' + repa.tagalign.basename + '_vs_' + \
                   repb.tagalign.basename
      proc = mp.Process(target=getNPeaksByIDR,
                        args=(fpeaka, fpeakb, idr_prefix, prm, out_q))
      procs.append(proc)
      proc.start()

    fidr2npeaks = {}
    for p in procs:
      fidr2npeaks.update(out_q.get())
      p.join()

    max_npeaks = max(fidr2npeaks.values())
    if not os.path.exists(self.peaks.fullname):
      sys.exit("File not found: %s\n" % self.peaks.fullname)

    with gzip.open(self.peaks.fullname, 'rb') as f_fin:
      sig_line = [ (float(line.split("\t")[6]), line) for line in f_fin ]
    sorted_sig_line = sorted(sig_line, key=lambda t: t[0], reverse=True)

    with gzip.open(self.final_peaks.fullname, 'wb') as f_fout:
      for (sig, line) in sorted_sig_line[:max_npeaks]:
        f_fout.write(line)


def getNPeaksByIDR(fpeaka, fpeakb, idr_prefix, prm, out_q):
  Util.runCommand('/bin/env', 'Rscript', prm.idr_script, fpeaka, fpeakb,
                  '-1', idr_prefix, '0', 'F', 'signal.value', prm.idr_scr_dir,
                  prm.fgenome_table, quiet=prm.quiet)
  fidr = idr_prefix + '-overlapped-peaks.txt'
  outdict = {}
  with open(fidr, 'r') as f_fidr:
    next(f_fidr)
    ## count the number of peaks w/ IDR <= IDR_THRESHOLD
    npk = sum( float(line.split()[10]) <= prm.IDR_THRESHOLD for line in f_fidr )
    outdict[fidr] = npk
    out_q.put(outdict)


def runSPP(tgt_tagalign, fctrl_tagalign, prm, nthr):
  spp_tmpdir = prm.temp_dir + tgt_tagalign.basename + '_spp_tmp/'
  if not os.path.exists(spp_tmpdir):
    os.mkdir(spp_tmpdir)
  fout = prm.temp_dir + tgt_tagalign.basename + '_phantom.tab'
  Util.runCommand('/bin/env', 'Rscript', prm.spp_script,
                  "-c=%s"      % tgt_tagalign.fullname,
                  "-i=%s"      % fctrl_tagalign,
                  "-npeak=%d"  % prm.N_PEAK,
                  prm.PEAK_TYPE,
                  '-savp',
                  "-x=%s"      % prm.EXCLUSION_ZONE,
                  '-rf',
                  "-odir=%s"   % prm.temp_dir,
                  "-p=%d"      % nthr,
                  "-tmpdir=%s" % spp_tmpdir,
                  "-out=%s"    % fout,
                  quiet=prm.quiet)
  Util.runCommand('rm', '-fr', spp_tmpdir, quiet=prm.quiet)

  if not os.path.exists(fout):
    sys.exit("Failed to generate file: %s\n" % fout)


def initFromParam(param, typ):
  if typ.lower() == 'target':
    is_ctrl = False
    param_attr = 'chipseq_target_read_files'
  elif typ.lower() in [ 'control', 'input' ]:
    is_ctrl = True
    param_attr = 'chipseq_control_read_files'
  elif typ.lower() == 'multi-targets':
    is_ctrl = False
    param_attr = 'chipseq_read_files_multi_targets'

  return ChIPSeqExperiment.initFromParam(param, is_ctrl, param_attr)
