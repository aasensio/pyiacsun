__all__ = ['tvframe']

import matplotlib.pyplot as pl
from mpl_toolkits.axes_grid1 import make_axes_locatable

def tvframe(image, ax=None, bar=True, **kwargs):
	"""Plot an image with a bar
	
	Args:
	    image (float): image
	    ax (axis): axis instance where to plot the image
	    bar (bool, optional): detault is to plot a bar. If not, set it to False
	    **kwargs: remaining arguments for the imshow
	
	Returns:
	    TYPE: None
	"""
	if (ax == None):
		f, ax = pl.subplots()
		
	im = ax.imshow(image, origin='lower', **kwargs)
	if (bar):
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.05)
		pl.colorbar(im, cax=cax)
