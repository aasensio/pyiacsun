
# Add fts_path to PATH
from . import ftsread 
import sys
ftsdir = str(ftsread.__file__).split('/')
strftsdir = ('/'.join(ftsdir[0:-1]))
sys.path.insert(0, strftsdir)
