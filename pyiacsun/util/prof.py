__all__ = ['prof4']
import pyfits as pf
import numpy as np
import matplotlib.pyplot as pl
import sys
import scipy.io

class prof4(object):
	def __init__(self, fileIn):
		self.currentFile = fileIn
		self.readData()
		self.currentPos = [0,0]
		
# Plot with integrated maps
		self.figFixedMaps = pl.figure(num=1, figsize=(5,8))
		self.axFixedMaps = [None]*4
		self.drawnFixedMaps = [None]*4
		self.slitMaps = self.getMaps(0)
		
		for i in range(4):
			self.axFixedMaps[i] = self.figFixedMaps.add_subplot(2,2,i+1)
			self.drawnFixedMaps[i] = self.axFixedMaps[i].imshow(self.maps[i], aspect='equal')
			self.axFixedMaps[i].set_axis_off()
			
		cid = self.figFixedMaps.canvas.mpl_connect('motion_notify_event', self.onMouseMove)
		cid = self.figFixedMaps.canvas.mpl_connect('close_event', self.onKillWindows)

# Plot with Stokes maps
		self.figMaps = pl.figure(num=2, figsize=(8,15))
		self.canvasMaps= self.figMaps.canvas
		self.axMaps = [None]*4
		self.drawnMaps = [None]*4
		self.posMaps = [None]*4
		self.backgroundMaps = [None]*4
		
		for i in range(4):
			self.axMaps[i] = self.figMaps.add_subplot(4,1,i+1)
			self.drawnMaps[i] = self.axMaps[i].imshow(self.slitMaps[i,:,:], aspect='equal', animated=True)			
			if (i == 0):
				self.posMaps[i], = self.axMaps[i].plot(200, 200, '+', animated=True)
			self.axMaps[i].set_xlim([0,self.nWavelength])
			self.axMaps[i].set_ylim([0,self.nx])
			
		self.canvasMaps.draw()
		
		for i in range(4):			
			self.backgroundMaps[i] = self.figMaps.canvas.copy_from_bbox(self.axMaps[i].bbox)
			
		cid = self.figMaps.canvas.mpl_connect('close_event', self.onKillWindows)
								
# Plot with Stokes profiles
		self.figStokes = pl.figure(num=3)
		self.canvasStokes = self.figStokes.canvas
		self.axStokes = [None]*4
		self.drawnStokes = [None]*4
		self.backgroundStokes = [None]*4
		cid = self.figStokes.canvas.mpl_connect('close_event', self.onKillWindows)
		
		for i in range(4):
			self.axStokes[i] = self.figStokes.add_subplot(2,2,i+1)
			self.drawnStokes[i], = self.axStokes[i].plot(self.slitMaps[i,0,:], animated=True)
		
		self.figStokes.tight_layout()
		self.canvasStokes.draw()
		
		self.newAxes = self.figStokes.get_axes()
		
		for i in range(4):
			self.axStokes[i] = self.newAxes[i]
			self.backgroundStokes[i] = self.figStokes.canvas.copy_from_bbox(self.axStokes[i].bbox)
			
	def updateProfiles(self, redrawMaps=False):
		for i in range(4):
			self.figStokes.canvas.restore_region(self.backgroundStokes[i])
			self.drawnStokes[i].set_ydata(self.slitMaps[i,self.currentPos[1],:])
			self.axStokes[i].draw_artist(self.drawnStokes[i])
			self.figStokes.canvas.blit(self.axStokes[i].bbox)
			
		if (redrawMaps):
			for i in range(4):
				self.figMaps.canvas.restore_region(self.backgroundMaps[i])
				self.drawnMaps[i].set_array(self.slitMaps[i,:,:])
				self.axMaps[i].draw_artist(self.drawnMaps[i])
				if (i == 0):
					self.posMaps[i].set_ydata(self.nx-self.currentPos[1])
					self.axMaps[i].draw_artist(self.posMaps[i])
					
				self.figMaps.canvas.blit(self.axMaps[i].bbox)
			
	def onMouseMove(self, event):	
		if (event.xdata != None and event.ydata != None):			
			newPos = np.asarray([event.xdata, event.ydata])
			newPos = newPos.astype(int)

			if (newPos[0] != self.currentPos[0]):
				self.slitMaps = self.getMaps(newPos[0])
				self.currentPos = newPos
				self.updateProfiles(redrawMaps=True)
			
			if (newPos[1] != self.currentPos[1]):
				self.currentPos = newPos
				self.updateProfiles(redrawMaps=False)
				
	def onKillWindows(self, event):
		pl.close('all')
		
	def readData(self):
				
		self.hdu = pf.open(self.currentFile+'c')
		self.nx = self.hdu[0].header['NAXIS2']
		self.nWavelength = self.hdu[0].header['NAXIS1']
		self.nFrames = self.hdu[0].header['NAXIS3'] / 4
				
		f = scipy.io.readsav(self.currentFile+'m')
		
		self.maps = [None] * 4
		nameMap = ['toti','totq','totu','totv']
		for i in range(4):
			self.maps[i] = np.flipud(f[nameMap[i]] / np.max(f[nameMap[i]]))
			
	
	def getProfiles(self, x, y):
		return self.hdu[0].data[4*x:4*x+4,y,:]
	
	def getMaps(self, x):
		return self.hdu[0].data[4*x:4*x+4,:,:]

#p = prof4(sys.argv[1])

#def main():
	#prof = prof4(sys.argv[1])
	
#if __name__ == "__main__":
   #main()
