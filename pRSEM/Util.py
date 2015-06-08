__doc__="""

  pliu 20150605

  utility module for pRSEM
  no class is defined here
"""

def runCommand(*args, **kwargs):
  import subprocess
  import sys

  str_args = [ str(arg) for arg in args ]
  if 'quiet' in kwargs:
    if not kwargs['quiet']:
      print ' '.join(str_args), "\n";
  else:
    print ' '.join(str_args), "\n";

  try:
    retcode = subprocess.call(str_args)
    if retcode < 0:
      sys.exit("Terminated by singal %d" % -retcode)
    elif retcode > 0:
      sys.exit("failed with return code %d" % retcode)
  except OSError as e:
    sys.exit("Execution failed: %s" % e)


def runCommandAndGetOutput(*args, **kwargs):
  import subprocess
  import sys

  str_args = [ str(arg) for arg in args ]
  if 'quiet' in kwargs:
    if not kwargs['quiet']:
      print ' '.join(str_args), "\n";
  else:
    print ' '.join(str_args), "\n";

  try:
    output = subprocess.check_output(str_args)
  except subprocess.CalledProcessError, e:
    sys.exit("Execution failed: %s" % e.output)

  return output


def runOneLineCommand(cmd, quiet=True):
  import os
  import sys

  if not quiet:
    print cmd, "\n";

  try:
    retcode = os.system(cmd)
    if retcode != 0:
      sys.exit("Failed with return code %d" % retcode)
  except OSError as e:
    sys.exit("Execution failed: %s" % e)


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
                         fbigwig):
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
