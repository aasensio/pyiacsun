import sys
from PyQt4 import QtCore, QtGui, uic
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
from ipdb import set_trace as stop
from netCDF4 import Dataset as nf
import matplotlib.pyplot as pl

 
form_class = uic.loadUiType("hazel.ui")[0]                 # Load the UI
 
class MyWindowClass(QtGui.QMainWindow, form_class):
	def __init__(self, fileIn, parent=None):
		QtGui.QMainWindow.__init__(self, parent)
		self.setupUi(self)

		self.currentFile = fileIn

		self.dpi = 80
		self.posMouse = [0,0]

		self.resize(1500,680)

		self.readData()
	
		self.leftFig = Figure((8,8), dpi=self.dpi)		
		self.leftCanvas = FigureCanvas(self.leftFig)
		self.leftCanvas.setParent(self.leftPlot)

		self.leftAxes = [None]*4
		self.drawnMap = [None]*4		
		for i in range(4):
			self.leftAxes[i] = self.leftFig.add_subplot(2,2,i+1)
			self.drawnMap[i] = self.leftAxes[i].imshow(self.maps[i], aspect='equal')

		self.leftFig.tight_layout()

		self.leftCanvas.mpl_connect('motion_notify_event', self.onMouseMove)
		self.leftCanvas.mpl_connect('scroll_event', self.onScrollMove)


		self.rightFig = Figure((8,8), dpi=self.dpi)		
		self.rightCanvas = FigureCanvas(self.rightFig)
		self.rightCanvas.setParent(self.rightPlot)

		self.stokesObs, self.stokesSyn, self.parsInv = self.getProfiles(0,0)
		self.axes = [None] * 4		
		self.obsLine = [None] * 4
		self.synLine = [None] * 4
		
		for i in range(4):
			self.axes[i] = self.rightFig.add_subplot(2, 2, i+1)
			self.obsLine[i], = self.axes[i].plot(self.stokesObs[:,i], animated=True)
			self.synLine[i], = self.axes[i].plot(self.stokesSyn[:,i], animated=True)
			self.axes[i].yaxis.set_animated(True)

		self.rightFig.tight_layout()
		
		self.rightCanvas.draw()
				
## We are using blit animation to do it fast, so we need to recover the axes modified by tight_layout
		self.newAxes = self.rightFig.get_axes()
		
## Save the backgrounds
		self.background = [None] * 4
		self.bbox = [None] * 4
		for i in range(4):
			self.axes[i] = self.newAxes[i]

			# self.bbox[i] = self.axes[i].bbox.union([label.get_window_extent() for label in self.axes[i].get_yticklabels()])			
			self.bbox[i] = self.axes[i].bbox.expanded(1.5, 1.1)
			self.background[i] = self.rightCanvas.copy_from_bbox(self.bbox[i])		

		self.tableWidget.setRowCount(9)
		self.tableWidget.setColumnCount(2)
		for i in range(9):
			self.tableWidget.setItem(i,0,QtGui.QTableWidgetItem("item {0}".format(i)))
			
 	def updateProfiles(self):		

# Blit animation. We only redraw the lines		
		
		for i in range(4):
			self.rightCanvas.restore_region(self.background[i])

			valMax = np.max([self.profilesObs[:,i],self.profilesSyn[:,i]])
			valMin = np.min([self.profilesObs[:,i],self.profilesSyn[:,i]])

# Redraw the observations
			self.obsLine[i].set_ydata(self.profilesObs[:,i])
			self.axes[i].draw_artist(self.obsLine[i])			
			
# Redraw the synthetic profiles
			self.synLine[i].set_ydata(self.profilesSyn[:,i])
			self.axes[i].draw_artist(self.synLine[i])

# Redraw the axis labels
			self.axes[i].set_ylim([valMin, valMax])
			self.axes[i].draw_artist(self.axes[i].get_yaxis())
		
			self.rightCanvas.blit(self.bbox[i])

			for i in range(9):
				item = QtGui.QTableWidgetItem("{0}".format(self.parsInv[i]))
				self.tableWidget.setItem(i, 1, item)

 	def onMouseMove(self, event):	
		if (event.xdata != None and event.ydata != None):			
			newPos = np.asarray([event.xdata, event.ydata])
			newPos = newPos.astype(int)

			if (newPos[0] != self.posMouse[0] or newPos[1] != self.posMouse[1]):
				self.posMouse = newPos
				self.profilesObs, self.profilesSyn, self.parsInv = self.getProfiles(newPos[0], newPos[1])
				self.updateProfiles()
			
	def onScrollMove(self,event):		
		pass

 	def getProfiles(self, x, y):
		return self.obs[y,x,:,0:5], self.syn[y,x,:,:], self.pars[y,x,:]
	
	def getMaps(self, x):
		return self.obs[:,:,0,0]

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


if __name__ == "__main__": 
	app = QtGui.QApplication(sys.argv)
	myWindow = MyWindowClass(sys.argv[1])
	myWindow.show()
	app.exec_()