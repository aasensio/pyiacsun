__all__ = ['display']
from Tkinter import Tk, Frame, Canvas
from PIL import Image, ImageTk
import numpy as np

class display(object):
	def __init__(self,image):
		self.data = np.array(image, dtype=int)
	
		self.root = Tk()
		self.frame = Tkinter.Frame(self.root, width=self.data.shape[1], height=self.data.shape[0])
		self.frame.pack()
		self.canvas = Tkinter.Canvas(self.frame, width=self.data.shape[1], height=self.data.shape[0])
		self.canvas.place(x=-2,y=-2)		
		self.im = Image.fromstring('L', (self.data.shape[1], self.data.shape[0]), self.data.astype('b').tostring())
		self.photo = ImageTk.PhotoImage(image=self.im)
		self.canvas.create_image(0,0,image=self.photo,anchor=Tkinter.NW)
		self.root.update()		
	
#im = 300*np.random.random((400,400))
#res = display(im)
	
