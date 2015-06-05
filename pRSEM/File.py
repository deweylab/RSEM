__doc__="""

  pliu 20150511

  modele for file-related definition and functions
"""

class File:
  def __init__(self):
    self.fullname = None  ## file's full name, include dir, base, and all ext
    self.is_gz    = None  ## if file is gzipped
    self.dirname  = None  ## directory name
    self.basename = None  ## base name sans all extension seperated by dot


  def __str__(self):
    ss = [ "fullname: %s\n" % self.fullname ] + \
         [ "dirname:  %s\n" % self.dirname  ] + \
         [ "basename: %s\n" % self.basename ]
    if self.is_gz:
      ss += [ "is gzipped" ]
    else:
      ss += [ "not gzipped" ]
    return ''.join(ss)


  @classmethod
  def initFromFullFileName(cls, filename):
    import os
    f = cls()
    f.fullname = filename
    (f.dirname, fname) = os.path.split(filename)
    words = fname.split('.')
    f.basename = words[0]
    if words[-1] in ['gz', 'gzip']:
      f.is_gz = True
    else:
      f.is_gz = False
    return f


def initFromFullFileName(ffq):
  return File.initFromFullFileName(ffq)
