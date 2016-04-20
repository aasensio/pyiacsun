from __future__ import print_function    
import pyiacsun

def test_atlas():
    atlas, xlam = pyiacsun.atlas.ftsread(ini = 10825, endi = 10843)
    atlas2, xlam2 = pyiacsun.atlas.Delbouille73(ini=6300, endi=6303)
    atlas3, xlam3 = pyiacsun.atlas.Neckel84(ini=6300, endi=6303)