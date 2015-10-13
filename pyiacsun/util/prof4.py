from __future__ import print_function

__all__ = ["AppForm"]

from PyQt4.QtGui import QMainWindow, QWidget, QApplication, QGridLayout
import sys
import os.path
	
import matplotlib.cm as cm
import numpy as np
import pyfits as pf
import scipy.io

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
class AppForm(QMainWindow):

	def __init__(self, fileIn, parent=None):
		QMainWindow.__init__(self, parent)
		
		self.setWindowTitle('PyProf')
		self.currentFile = fileIn
				
		self.whichShowing = 0
		self.posMouse = [0,0]

		self.create_main_frame()
				
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
		print(event.xdata, event.ydata)
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
			
	def readData(self):
				
		self.hdu = pf.open(self.currentFile+'c')
		self.nx = self.hdu[0].header['NAXIS2']
		self.nWavelength = self.hdu[0].header['NAXIS1']
		self.nFrames = self.hdu[0].header['NAXIS3'] / 4
				
		f = scipy.io.readsav(self.currentFile+'m')
		
		self.maps = [None] * 4
		nameMap = ['toti','totq','totu','totv']
		for i in range(4):
			self.maps[i] = f[nameMap[i]] / np.max(f[nameMap[i]])			
	
	def getProfiles(self, x, y):
		return self.hdu[0].data[4*x:4*x+4,y,:]
	
	def getMaps(self, x):
		return self.hdu[0].data[4*x:4*x+4,:,:]
		
		
#--------------------------------------    
	def create_main_frame(self):
        
		self.mainFrame = QWidget()
		
		gridLayout = QGridLayout()
                       
		self.dpi = 80
		self.xpos = 0
		self.ypos = 0
				
		self.titles = [r'I',r'Q',r'U',r'V']
		
		self.readData()
		
		self.resize(1500,500)
			
# Left window
		self.leftPlot = QWidget()
		#self.leftFig = Figure(((4*self.nFrames) / self.dpi,2*self.nx / self.dpi), dpi=self.dpi)
		self.leftFig = Figure(((4*self.nFrames) / self.dpi,8), dpi=self.dpi)
		self.leftCanvas = FigureCanvas(self.leftFig)
		self.leftCanvas.setParent(self.leftPlot)
					
		self.leftAxes = [None]*4
		self.drawnMap = [None]*4
		for i in range(4):
			self.leftAxes[i] = self.leftFig.add_subplot(2,2,i+1)
			self.drawnMap[i] = self.leftAxes[i].imshow(self.maps[i], aspect='equal')
			self.leftAxes[i].set_axis_off()
		self.leftCanvas.mpl_connect('motion_notify_event', self.onMouseMove)
		self.leftCanvas.mpl_connect('scroll_event', self.onScrollMove)		
		gridLayout.addWidget(self.leftPlot, 0, 0)
		gridLayout.setSpacing(10)
		
# Central window				
		self.centralPlot = QWidget()		
		#self.centralFig = Figure((self.nWavelength / self.dpi,2*self.nx / self.dpi), dpi=self.dpi)
		self.centralFig = Figure((10,8), dpi=self.dpi)
		self.centralCanvas = FigureCanvas(self.centralFig)
		self.centralCanvas.setParent(self.centralPlot)
				
		self.slitMaps = self.getMaps(0)
		
		self.centralAxes = [None]*4
		self.drawnSlitMap = [None]*4
		for i in range(4):
			self.centralAxes[i] = self.centralFig.add_subplot(4,1,i+1)
			self.drawnSlitMap[i] = self.centralAxes[i].imshow(self.slitMaps[i,:,:])
		#self.centralCanvas.draw()
			
		gridLayout.addWidget(self.centralPlot, 0, 1)

# Right window		
		self.rightPlot = QWidget()
		self.rightFig = Figure((8,8), dpi=self.dpi)
		self.rightCanvas = FigureCanvas(self.rightFig)
		self.rightCanvas.setParent(self.rightPlot)
							
		self.stokes = self.getProfiles(0,0)
		
# Draw the axes and the first profiles
		nCols = 2
		nRows = 2
		self.axes = [None] * 4
		self.obs = [None] * 4		
		loop = 0
		for j in range(2):
			for i in range(2):
				self.axes[loop] = self.rightFig.add_subplot(nCols, nRows, loop+1)
				self.obs[loop], = self.axes[loop].plot(self.stokes[loop,:])
				loop += 1
		gridLayout.addWidget(self.rightPlot, 0, 2)		
		
		
## Tight layout and redraw
		self.rightFig.tight_layout()
		
		self.rightCanvas.draw()
				
## We are using blit animation to do it fast, so we need to recover the axes modified by tight_layout
		self.newAxes = self.figRight.get_axes()
		
## Save the backgrounds
		loop = 0
		for j in range(6):
			for i in range(4):
				self.axes[loop] = self.newAxes[loop]
				self.background[loop] = self.canvasRight.copy_from_bbox(self.axes[loop].bbox)
				loop += 1
															
		
		
		
		plotLayout.addWidget(self.rightPlot)
		
		fullLayout.addLayout(gridLayout)
		
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