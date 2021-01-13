# imports
import os
import time
from threading import Thread
from whitebox.ui.plugin_dialog import ScriptDialog
from java.awt.event import ActionListener
from whitebox.geospatialfiles import WhiteboxRaster
from whitebox.geospatialfiles.WhiteboxRasterBase import DataType

'''The following four variables are required for this 
   script to be integrated into the tool tree panel. 
   Comment them out if you want to remove the script.'''
#name = "PythonExamplePlugin" 
descriptiveName = "Example Python Plugin" 
description = "Just an example of a plugin tool using Python." 
toolboxes = ["topmost"] 
	
class PythonExamplePlugin(ActionListener):
	def __init__(self, args):
		if len(args) != 0:
			self.execute(args)
		else:
			''' Create a dialog for this tool to collect user-specified
			   tool parameters.''' 
			self.sd = ScriptDialog(pluginHost, "Python Example Plugin", self)	
			
			''' Specifying the help file will display the html help
			// file in the help pane. This file should be be located 
			// in the help directory and have the same name as the 
			// class, with an html extension.'''
			helpFile = self.__class__.__name__
			self.sd.setHelpFile(helpFile)
	
			''' Specifying the source file allows the 'view code' 
			// button on the tool dialog to be displayed.'''
			self.sd.setSourceFile(os.path.abspath(__file__))
	
			# add some components to the dialog '''
			self.sd.addDialogFile("Input raster file", "Input Raster File:", "open", "Raster Files (*.dep), DEP", True, False)
			self.sd.addDialogFile("Output raster file", "Output Raster File:", "save", "Raster Files (*.dep), DEP", True, False)
			
			# Resize the dialog to the standard size and display it '''
			self.sd.setSize(800, 400)
			self.sd.visible = True

	def actionPerformed(self, event):
		if event.getActionCommand() == "ok":
			args = self.sd.collectParameters()
			t = Thread(target=lambda: self.execute(args))
			t.start()

	''' The execute function is the main part of the tool, where the actual
        work is completed.'''
	def execute(self, args):
		try:
			dX = [ 1, 1, 1, 0, -1, -1, -1, 0 ]
			dY = [ -1, 0, 1, 1, 1, 0, -1, -1 ]
			
			if len(args) != 2:
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return

			# read the input parameters
			inputfile = args[0]
			outputfile = args[1]
			
			# read the input image 
			inputraster = WhiteboxRaster(inputfile, 'r')
			nodata = inputraster.getNoDataValue()
			rows = inputraster.getNumberRows()
			cols = inputraster.getNumberColumns()
			
			# initialize the output image
			outputraster = WhiteboxRaster(outputfile, "rw", inputfile, DataType.FLOAT, nodata)
			outputraster.setPreferredPalette(inputraster.getPreferredPalette())

			'''perform the analysis
			   This code loops through a raster and performs a 
		 	   3 x 3 mean filter.'''
			oldprogress = -1
			for row in xrange(0, rows):
				for col in xrange(0, cols):
					z = inputraster.getValue(row, col)
					if z != nodata:
						mean = z
						numneighbours = 1
						for n in xrange(0, 8):
							zn = inputraster.getValue(row + dY[n], col + dX[n])
							if zn != nodata:
								mean += zn
								numneighbours += 1
								
						outputraster.setValue(row, col, mean / numneighbours)
				
				progress = (int)(100.0 * row / (rows - 1))
				if progress > oldprogress:
					oldprogress = progress
					pluginHost.updateProgress(progress)
					if pluginHost.isRequestForOperationCancelSet():
						pluginHost.showFeedback("Operation cancelled")
						return
			
			inputraster.close()
			outputraster.addMetadataEntry("Created by the " + descriptiveName + " tool.")
			outputraster.addMetadataEntry("Created on " + time.asctime())
			outputraster.close()

			# display the output image
			pluginHost.returnData(outputfile)
			
		except Exception, e:
			print e
			pluginHost.showFeedback("An error has occurred during operation. See log file for details.")
			pluginHost.logException("Error in " + descriptiveName, e)
			return
		finally:
			# reset the progress bar
			pluginHost.updateProgress(0)
			
if args is None:
	pluginHost.showFeedback("The arguments array has not been set.")
else:		
	PythonExamplePlugin(args)
