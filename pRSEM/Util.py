__doc__="""

  pliu 20150605

  utility module for pRSEM
  no class is defined here
"""

def runCommand(*args):
  import subprocess
  import sys

  str_args = [ str(arg) for arg in args ]
  print ' '.join(str_args), "\n";
  try:
    retcode = subprocess.call(str_args)
    if retcode < 0:
      sys.exit("Terminated by singal %d" % -retcode)
    elif retcode > 0:
      sys.exit("failed with return code %d" % retcode)
  except OSError as e:
    sys.exit("Execution failed: %s" % e)


def runOneLineCommand(cmd):
  import os
  import sys

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
