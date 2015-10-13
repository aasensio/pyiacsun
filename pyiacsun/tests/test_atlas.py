from __future__ import print_function    
import pyiacsun

def test_atlas():
    atlas, xlam = pyiacsun.atlas.ftsread(ini = 10825, endi = 10843)