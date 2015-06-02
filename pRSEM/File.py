__doc__="""

  pliu 20150511

  modele for file-related definition and functions
"""

class File:
  def __init__(self):
    self.fullname = None  ## file's full name, include dir, base, and all ext
    self.is_gz    = None  ## if file is gzipped
    self.dirname  = None  ## directory name
    self.basename = None  ## base name san all extension seperated by .
    self.type     = None

  @classmethod
  def initFromFullFileName(cls, ffq):
    f = cls()
    return f


def initFromFullFileName(ffq):
  return File.initFromFullFileName(ffq)
