__all__ = ['cursor']

from matplotlib.widgets import Cursor
import matplotlib.pyplot as pl
import numpy as np

class cursor(object):
	def __init__(self, fig):
		self.fig = fig
		self.ax = self.fig.get_axes()		
		c = Cursor(self.ax[0], useblit=True)
		self._widgets=[c]
		
		self.x = []
		self.y = []
		
		self.connectID = self.fig.canvas.mpl_connect('button_press_event', self.click)
				
# start a blocking event loop
		self.fig.canvas.start_event_loop(timeout=-1)
		
	def click(self, event):
		if (event.button == 1):
			if (event.inaxes == self.ax[0]):				
				print(event.xdata, event.ydata)
				self.x.append(event.xdata)
				self.y.append(event.ydata)
		if (event.button == 3):
			self.fig.canvas.mpl_disconnect(self.connectID)
			self.fig.canvas.stop_event_loop()
		
		
#f = pl.figure(num=1)
#ax = f.add_subplot(1,1,1)
#ax.plot(np.linspace(1,10,100))

#res = cursor(f)