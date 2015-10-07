# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 10:18:41 2014

@author: aasensio
"""
__all__ = ["AppForm"]

from PyQt4.QtGui import QMainWindow, QWidget, QApplication, QGridLayout
import sys
	
import numpy as np
from netCDF4 import Dataset as nf
from ipdb import set_trace as stop

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

class AppForm(QMainWindow):

	def __init__(self, fileIn, parent=None):
		QMainWindow.__init__(self, parent)
		
		self.setWindowTitle('PyProf')
		self.currentFile = fileIn
				
		self.whichShowing = 0
		self.posMouse = [0,0]

		self.createMainFrame()

	def readParameters(self, rootFile):
	
# Read inverted parameters
		ff = nf(rootFile+'.parameters', 'r')
		pars = ff.variables['map'][:]
		ff.close()
	
# Read errors
		ff = nf(rootFile+'.errors', 'r')
		errors = ff.variables['map'][:]
		ff.close()
		
		return pars, errors

	def readProfiles(self, rootFile):

# Read inverted profiles
		ff = nf(rootFile+'.inversion', 'r')
		synthProf = ff.variables['map'][:]
		ff.close()
	
		return synthProf

	def readObservations(self, rootFile):

# Read inverted profiles
		ff = nf(rootFile+'.nc', 'r')
		sizeMask = ff.variables['mask'].shape
		sizeMap = ff.variables['map'].shape
		obsProf = ff.variables['map'][:].reshape((sizeMask[0],sizeMask[1],sizeMap[-2],sizeMap[-1]))
		ff.close()
	
		return obsProf
			
	def readData(self):
	
	# Parameters
		self.obs = self.readObservations(self.currentFile)
		obsShape = self.obs.shape

		self.pars, self.errors = self.readParameters(self.currentFile)
		parsShape = self.pars.shape
		self.pars = self.pars.reshape((obsShape[0],obsShape[1],parsShape[-1]))
		self.errors = self.errors.reshape((obsShape[0],obsShape[1],parsShape[-1]))
		
		self.syn = self.readProfiles(self.currentFile)
		synShape = self.syn.shape
		self.syn = self.syn.reshape((obsShape[0],obsShape[1],synShape[-2],synShape[-1]))
				
		self.nx, self.ny, self.nLambda, _ = self.syn.shape		
						
		self.maps = [None] * 4
		self.maps = [None] * 4		
		for i in range(4):
			self.maps[i] = np.sum(self.obs[:,:,:,0], axis=(2))
				
	def updateProfiles(self):		
# Blit animation. We only redraw the lines
		loop = 0
		for j in range(6):
			for i in range(4):
				self.rightCanvas.restore_region(self.background[loop])
				self.obs[loop].set_ydata(self.profiles[self.rangeFrom[j]:self.rangeTo[j],i])
				self.axes[loop].draw_artist(self.obs[loop])
				self.rightCanvas.blit(self.axes[loop].bbox)
				loop += 1		
						
	def onMouseMove(self, event):	
		print event.xdata, event.ydata
		if (event.xdata != None and event.ydata != None):			
			newPos = np.asarray([event.xdata, event.ydata])
			newPos = newPos.astype(int)

			if (newPos[0] != self.posMouse[0] or newPos[1] != self.posMouse[1]):
				self.posMouse = newPos
				self.profiles = self.getProfiles(newPos[0], newPos[1])
				self.updateProfiles()
			
	def onScrollMove(self,event):		
		pass
		#if (event.button == 'up'):
			#self.newTime += 1			
		#if (event.button == 'down'):
			#self.newTime -= 1
			
		#self.newTime = self.newTime % 35
		
		#self.profiles = self.getProfiles(self.posMouse[0], self.posMouse[1])
		#self.updateProfiles()

			
	def getProfiles(self, x, y):
		return self.obs[x,y,:,0:5], self.syn[x,y,:,:]
	
	def getMaps(self, x):
		return self.obs[:,:,0,0]
		
		
#--------------------------------------    
	def create_main_frame(self):
        
		self.mainFrame = QWidget()
		
		gridLayout = QGridLayout()
                       
		self.dpi = 80
		self.xpos = 0
		self.ypos = 0
				
		self.titles = [r'I',r'Q',r'U',r'V']
		
		self.readData()
		
		self.resize(1500,900)
			
# Left window
		self.leftPlot = QWidget()		
		self.leftFig = Figure((2*self.ny / self.dpi, 4*self.nx / self.dpi), dpi=self.dpi)		
		self.leftCanvas = FigureCanvas(self.leftFig)
		self.leftCanvas.setParent(self.leftPlot)
					
		self.leftAxes = [None]*4
		self.drawnMap = [None]*4		
		for i in range(4):
			self.leftAxes[i] = self.leftFig.add_subplot(2,2,i+1)
			self.drawnMap[i] = self.leftAxes[i].imshow(self.maps[i], aspect='equal')
			# self.leftAxes[i].set_axis_off()
		self.leftCanvas.mpl_connect('motion_notify_event', self.onMouseMove)
		self.leftCanvas.mpl_connect('scroll_event', self.onScrollMove)		
		gridLayout.addWidget(self.leftPlot, 0, 0)
		gridLayout.setSpacing(10)
		
# Central window				
		# self.centralPlot = QWidget()		
		# self.centralFig = Figure((10,8), dpi=self.dpi)
		# self.centralCanvas = FigureCanvas(self.centralFig)
		# self.centralCanvas.setParent(self.centralPlot)
				
		# self.slitMaps = self.getMaps(0)
		
		# self.centralAxes = [None]*4
		# self.drawnSlitMap = [None]*4
		# # for i in range(4):
		# # 	self.centralAxes[i] = self.centralFig.add_subplot(4,1,i+1)
		# # 	self.drawnSlitMap[i] = self.centralAxes[i].imshow(self.slitMaps[i,:,:])
		# #self.centralCanvas.draw()
			
		# gridLayout.addWidget(self.centralPlot, 0, 1)

# Right window		
		self.rightPlot = QWidget()
		self.rightFig = Figure((8,18), dpi=self.dpi)
		self.rightCanvas = FigureCanvas(self.rightFig)
		self.rightCanvas.setParent(self.rightPlot)
							
		self.stokesObs, self.stokesSyn = self.getProfiles(0,0)
		
# Draw the axes and the first profiles
		nCols = 2
		nRows = 2
		self.axes = [None] * 4
		self.obsLine = [None] * 4
		self.synLine = [None] * 4
		loop = 0
		for j in range(2):
			for i in range(2):
				self.axes[loop] = self.rightFig.add_subplot(nCols, nRows, loop+1)
				self.obsLine[loop], = self.axes[loop].plot(self.stokesObs[:,loop])
				self.synLine[loop], = self.axes[loop].plot(self.stokesSyn[:,loop])
				loop += 1
		gridLayout.addWidget(self.rightPlot, 0, 2)		
		
		
## Tight layout and redraw
		# self.rightFig.tight_layout()
		
		self.rightCanvas.draw()
				
## We are using blit animation to do it fast, so we need to recover the axes modified by tight_layout
		self.newAxes = self.rightFig.get_axes()
		
## Save the backgrounds
		# loop = 0
		# for j in range(6):
		# 	for i in range(4):
		# 		self.axes[loop] = self.newAxes[loop]
		# 		self.background[loop] = self.canvasRight.copy_from_bbox(self.axes[loop].bbox)
		# 		loop += 1
															
					
		gridLayout.addWidget(self.rightPlot)
		
		# fullLayout.addLayout(gridLayout)
		
		self.mainFrame.setLayout(gridLayout)
		self.setCentralWidget(self.mainFrame)
		
class SideForm(QMainWindow):
	def __init__(self, fileIn, parent=None):
		QMainWindow.__init__(self, parent)
		
def main():
	app = QApplication(sys.argv)
	form = AppForm(sys.argv[1])
	form.show()
	app.exec_()
	del form
	
if __name__ == "__main__":
   main()