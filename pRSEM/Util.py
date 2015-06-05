__doc__="""

  pliu 20150605

  utility module for pRSEM
  no class is defined here
"""

def runCommand(*args):
  import subprocess
  import sys

  str_args = [ str(arg) for arg in args ]
  cmd = ' '.join(str_args)
  print cmd;
  try:
    retcode = subprocess.call(str_args)
    if retcode < 0:
      print >> sys.stderr, 'Child was terminated by singal', -retcode;
  except OSError as e:
    print >> sys.stderr, 'Execution failed:', e;


def getCatCommand(is_gzipped):
  if is_gzipped:
    cat_cmd = 'zcat'
  else:
    cat_cmd = 'cat'
  return cat_cmd
