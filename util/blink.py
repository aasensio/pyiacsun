import numpy as np
import matplotlib.pyplot as pl
import matplotlib.animation as animation

def blink(ima1, ima2, interval=250, fig=None):
	"""Blink two images
	
	Args:
	    ima1 (float): first image
	    ima2 (float): second image
	    interval (int, optional): time step to change from one to the other
	    fig (TYPE, optional): optional figure handle to do the blinking in an existing figure
	
	Returns:
	    TYPE: None
	"""
	which = 0
	if (fig == None):
		fig = pl.figure()
	im = pl.imshow(ima1)
	print(which)

	def updatefig(*args):
		global which
		if (which == 0):
			im.set_array(args[0])
			which = 1
		else:
			im.set_array(args[1])
			which = 0
		return im,

	ani = animation.FuncAnimation(fig, updatefig, interval=interval, blit=True, fargs=(ima1,ima2))