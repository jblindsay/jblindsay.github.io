'''
 * Copyright (C) 2013 Dr. John Lindsay <jlindsay@uoguelph.ca>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
 
import os
from shutil import copyfile
from threading import Thread
from whitebox.geospatialfiles import WhiteboxRaster
from whitebox.ui.plugin_dialog import ScriptDialog
from java.awt.event import ActionListener
from java.awt.event import ActionEvent

# The following four variables are required for this 
# script to be integrated into the tool tree panel. 
# Comment them out if you want to remove the script.
name = "Sink"
descriptiveName = "Sink"
description = "Identifies the sinks (depressions) in a DEM"
toolboxes = ["TerrainAnalysis"]
	
class Sink(ActionListener):
	def __init__(self, args):
		if len(args) != 0:
			self.execute(args)
		else:
			# Create a dialog for this tool to collect user-specified
			# tool parameters.
			self.sd = ScriptDialog(pluginHost, descriptiveName, self)	
		
			# Specifying the help file will display the html help
			# file in the help pane. This file should be be located 
			# in the help directory and have the same name as the 
			# class, with an html extension.
			helpFile = self.__class__.__name__
			self.sd.setHelpFile(helpFile)
		
			# Specifying the source file allows the 'view code' 
			# button on the tool dialog to be displayed.
			self.sd.setSourceFile(os.path.abspath(__file__))

			# add some components to the dialog
			self.sd.addDialogFile("Input file", "Input DEM File:", "open", "Whitebox Raster Files (*.dep), DEP", True, False)
			self.sd.addDialogFile("Output file", "Output Raster File:", "close", "Whitebox Raster Files (*.dep), DEP", True, False)

			# resize the dialog to the standard size and display it
			self.sd.setSize(800, 400)
			self.sd.visible = True
		
	def execute(self, args):
		try:
			inputfile = args[0]
			outputfile = args[1]

			# fill the depressions in the DEM, being sure to 
			# run on a dedicated thread and supressing
			# the return data (automatically displayed image)
			outputfile2 = outputfile.replace(".dep", "_temp.dep")
			args2 = [inputfile, outputfile2, "0.0"]
			pluginHost.runPlugin("FillDepressions", args2, False, True)

			# flag the cells that have been affected by the filling.
			inputraster = WhiteboxRaster(inputfile, 'r')
			tempraster = WhiteboxRaster(outputfile2, 'rw')
			rows = tempraster.getNumberRows()
			cols = tempraster.getNumberColumns()
			nodata = inputraster.getNoDataValue()

			oldprogress = -1
			for row in xrange(0, rows):
				for col in xrange(0, cols):
					z1 = inputraster.getValue(row, col)
					z2 = tempraster.getValue(row, col)
					if z1 != nodata:
						if z1 < z2:
							tempraster.setValue(row, col, 1.0)
						else:
							tempraster.setValue(row, col, 0.0)
					else:
						tempraster.setValue(row, col, nodata)
				
				progress = (int)((100.0 * (row + 1.0)) / rows)
				if progress > oldprogress:
					oldprogress = progress
					pluginHost.updateProgress(progress)
				if pluginHost.isRequestForOperationCancelSet():
					pluginHost.showFeedback("Operation cancelled")
					return
					
			tempraster.flush()

			# clump the filled cells, being sure to 
			# run on a dedicated thread and supressing
			# the return data (automatically displayed image)
			args2 = [outputfile2, outputfile, 'true', 'true']
			pluginHost.runPlugin("Clump", args2, False, True)

			inputraster.close()
			
			# delete the temporary file
			os.remove(outputfile2)
			os.remove(outputfile2.replace(".dep", ".tas"))
			
			# display the final output
			pluginHost.returnData(outputfile)
			
		except:
			pluginHost.showFeedback("Error during script execution.")
			return

	def actionPerformed(self, event):
		if event.getActionCommand() == "ok":
			args = self.sd.collectParameters()
			self.sd.dispose()
			t = Thread(target=lambda: self.execute(args))
			t.start()

if args is None:
	pluginHost.showFeedback("The arguments array has not been set.")
else:	
	Sink(args)
