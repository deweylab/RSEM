__doc__="""

  pliu 20150605

  utility module for pRSEM
  no class is defined here
"""

def runCommand(*args, **kwargs):
  import os
  import subprocess
  import sys

  is_quiet = False
  if 'quiet' in kwargs:
    if kwargs['quiet']:
      is_quiet = True

  str_args = [ str(arg) for arg in args ]
  if is_quiet:
    pass
  else:
    sys.stdout.write("\n%s\n" % (' '.join(str_args)))

  f_null = open(os.devnull, 'w')

  try:
    if len(str_args) == 1:
      if is_quiet:
        retcode = subprocess.call(str_args[0], stdout=f_null, shell=True)
      else:
        retcode = subprocess.call(str_args[0], shell=True)
    else:
      if is_quiet:
       #print '##', is_quiet, '##';
        retcode = subprocess.call(str_args, stdout=f_null)
      else:
       #print '##', is_quiet, '##';
        retcode = subprocess.call(str_args)
    if retcode < 0:
      sys.exit("\nTerminated by singal %d\n" % -retcode)
    elif retcode > 0:
      sys.exit("\nFailed with return code %d\n" % retcode)
  except OSError as e:
    sys.exit("\nExecution failed: %s\n" % e)

  f_null.close()


def runCommandAndGetOutput(*args, **kwargs):
  import subprocess
  import sys

  str_args = [ str(arg) for arg in args ]
  if 'quiet' in kwargs:
    if not kwargs['quiet']:
      sys.stdout.write("\n%s\n" % (' '.join(str_args)))
  else:
    sys.stdout.write("\n%s\n" % (' '.join(str_args)))

  try:
    output = subprocess.check_output(str_args)
  except subprocess.CalledProcessError, e:
    sys.exit("\nExecution failed: %s\n" % e.output)

  return output


def getCatCommand(is_gzipped):
  if is_gzipped:
    cat_cmd = 'zcat'
  else:
    cat_cmd = 'cat'
  return cat_cmd


def readFile(fin):
  """
  return all the lines of the input file.
  """
  import os
  assert os.path.exists(fin), "File not found: %s\n" % fin

  lines = [];
  f_fin = open(fin, 'r');
  lines = f_fin.read().split('\n');
  f_fin.close();
  lines.pop();

  newlines = [];
  for line in lines:
    if line[-1] == '\r':
      newline = line[:-1];
    else:
      newline = line;
    newlines.append(newline);

  return newlines;


def calculateMappability(mpp_type, chrom, start, end, bigwigsummary_bin,
                         fbigwig, quiet=True):
  """
  calculate mappability for the given genomic coordinate interval
  mpp_type = {mean|max}
  """
  mpp = -10.0
  mpp = runCommandAndGetOutput(bigwigsummary_bin, '-type=%s' % mpp_type,
                               fbigwig, chrom, start, end, '1', quiet=True)
  return float(mpp)


def runMPOverAList(nprocs, func, args):
  """
  run multiprocessing for the given function and arguments on nprocs CPUs
  args[0] must be a list to-be-split and run func
  func must return a dict
  """
  import multiprocessing as mp

  out_q = mp.Queue()
  chunksize = 1
  if len(args[0]) > nprocs:
    chunksize = len(args[0])/nprocs + 1
  procs = []
  for i in xrange(nprocs):
    list_args = [args[0][chunksize*i:chunksize*(i+1)]] + args[1:] + [out_q]
    p = mp.Process(target = func, args = tuple(list_args))
    procs.append(p)
    p.start()

  dict_to_return = {}
  for i in xrange(nprocs):
    dict_to_return.update(out_q.get())

  for p in procs:
    p.join()

  return dict_to_return


def getFastaID2Seq(ffasta):
  """
  read fasta file, return a dict with key as seq_id and value as seq
  """
  import os
  assert os.path.exists(ffasta), "File not found: %s\n" % ffasta
  fastas = {};
  f_fin = open(ffasta, 'r');
  entries = f_fin.read().split('>');
  f_fin.close();
  for entry in entries[1:]:
    words = entry.split("\n");
    fastas[words[0]] = words[1];

  return fastas;


def getGCFraction(seq):
  """
  return the percetage of GC in the given sequence
  """
  length = len(seq);
  if length == 0:
    sys.stderr.write("Util::getFraction(): sequence length is 0\n");
    return 0;
  else:
    seq = seq.upper();
    n_G = seq.count('G');
    n_C = seq.count('C');

    return (n_G + n_C) * 1.0/length;

