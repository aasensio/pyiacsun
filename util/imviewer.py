#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Author: cdiazbas@iac.es
# Date: 28.04.2015
# Code: Gregor Images Viewer



from matplotlib.widgets import Cursor, Button
import matplotlib.pyplot as plt
import numpy as np
# plt.rcParams['font.size']=8

class imviewer(object):

    def __init__(self,imagen,stokeim=0,stokeplot=0):

        # Stoke in plot:
        globals()['stokeplot'] = stokeplot
        # Stoke in imshow:
        globals()['stokeim'] = stokeim

        print('Click derecho para borrar todo.')
        self.z=imagen[:,stokeim,:,0]
        lenTime=imagen.shape[0]
        lenStokes=imagen.shape[1]
        lenSlit=imagen.shape[2]
        lenLambda=imagen.shape[3]
        self.x=np.arange(lenSlit)
        self.y=np.arange(lenTime)

        self.fig=plt.figure(figsize=(8, 6))
        
        #Doing some layout with subplots:
        self.fig.subplots_adjust(0.08,0.08,0.98,0.98,0.1)
        self.overview=plt.subplot2grid((2,2),(0,0),rowspan=1,colspan=3)
        self.overview.imshow(self.z, cmap='gray',origin='lower',interpolation='None')
        plt.ylabel('Time Axis')
        plt.xlabel('Slit Axis')
        plt.ylim((0,lenTime-1))
        plt.xlim((0,lenSlit-1))

        self.overview.autoscale(1,'both',1)
        self.mini_subplot=plt.subplot2grid((2,2),(1,0),rowspan=1,colspan=3)
        plt.ylabel('Stokes value')
        plt.xlabel('Lambda')
        plt.grid(True)
        plt.ylim((0,1.2))
        plt.xlim((0,imagen.shape[3]-1))

        #Adding widgets, to not be gc'ed, they are put in a list:
        cursor=Cursor(self.overview, useblit=True, color='black', linewidth=1 )
        cursor2=Cursor(self.mini_subplot, useblit=True, color='black', linewidth=1 )
        self._widgets=[cursor,cursor2]
        
        #connect events
        self.fig.canvas.mpl_connect('button_press_event',self.click)
    


    def click(self,event):
        """
        What to do, if a click on the figure happens:
            1. Check which axis
            2. Get data coord's.
            3. Plot resulting data.
            4. Update Figure
        """

        if event.inaxes==self.overview:
            
            #Get nearest data
            xpos=np.argmin(np.abs(event.xdata-self.x))
            ypos=np.argmin(np.abs(event.ydata-self.y))
            nimagen = imagen[ypos,stokeplot,xpos,:]/1e+4
            
            #Check which mouse button:
            if event.button==1:
                print(xpos,ypos)      
                c,=self.mini_subplot.plot(np.arange(imagen.shape[3]),nimagen,label=str(self.x[xpos]))
                self.overview.scatter(self.x[xpos],self.y[ypos],color=c.get_color())
                self.overview.set_ylim((0,imagen.shape[0]-1))
                self.overview.set_xlim((0,imagen.shape[2]-1))

            elif event.button==3:
                
                # Clear the plot
                plt.clf()
                #Doing some layout with subplots:
                self.fig.subplots_adjust(0.08,0.08,0.98,0.98,0.1)
                self.overview=plt.subplot2grid((2,2),(0,0),rowspan=1,colspan=3)
                self.overview.imshow(self.z, cmap='gray',origin='lower')
                plt.ylabel('Time Axis')
                plt.xlabel('Slit Axis')
                plt.ylim((0,imagen.shape[0]-1))
                plt.xlim((0,imagen.shape[2]-1))

                self.overview.autoscale(1,'both',1)
                self.mini_subplot=plt.subplot2grid((2,2),(1,0),rowspan=1,colspan=3)
                plt.ylabel('Stokes value')
                plt.xlabel('Lambda')
                plt.grid(True)
                plt.ylim((0,1.2))
                plt.xlim((0,imagen.shape[3]-1))

                #Adding widgets, to not be gc'ed, they are put in a list:
                cursor=Cursor(self.overview, useblit=True, color='black', linewidth=1 )
                cursor2=Cursor(self.mini_subplot, useblit=True, color='black', linewidth=1 )
                self._widgets=[cursor,cursor2]

        if event.inaxes==self.mini_subplot:
            self.x=np.arange(imagen.shape[3])
            xpos=np.argmin(np.abs(event.xdata-self.x))
            nimagen = imagen[:,stokeim,:,xpos]
            c=self.overview.imshow(nimagen, cmap='gray',origin='lower',interpolation='None')
            plt.ylabel('Stokes value')
            plt.xlabel('Lambda')
            plt.grid(True)
            plt.ylim((0,1.2))
            plt.xlim((0,imagen.shape[3]-1))


        #Show it
        plt.draw()
        
if __name__=='__main__':

    import cjd_pylib.imtools as imtools

    FileName = '17jun14.006-01ccms.fits'

    imagen, header = imtools.imfits(FileName)
    # imagen, header = imtools.imgregor(FileName,info = False)

    fig_v=imviewer(imagen,stokeim=3,stokeplot=0)
    plt.show()


