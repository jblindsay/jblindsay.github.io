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
from whitebox.ui.plugin_dialog import ScriptDialog
from whitebox.geospatialfiles import GeoTiff
from whitebox.geospatialfiles import GeoKey
from whitebox.geospatialfiles import IFDEntry
from java.awt.event import ActionListener
from java.awt.event import ActionEvent

# The following four variables are required for this 
# script to be integrated into the tool tree panel. 
# Comment them out if you want to remove the script.
name = "ShowGeoTiffTags"
descriptiveName = "Show GeoTiff Tags"
description = "Prints the tags belonging to a GeoTiff"
toolboxes = ["FileUtilities", "IOTools"]
	
class ShowGeoTiffTags(ActionListener):
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
			self.sd.addDialogFile("Input file", "Input GeoTiff File:", "open", "GeoTiff Files (*.tif), TIF, TIFF", True, False)
			
			# resize the dialog to the standard size and display it
			self.sd.setSize(800, 400)
			self.sd.visible = True
		
	def execute(self, args):
		try:
			inputfile = args[0]

			if not os.path.exists(inputfile):
				pluginHost.showFeedback("File " + inputfile + " does not exist.")
				return
					
			gt = GeoTiff(inputfile)
			gt.read()
			tags = gt.showInfo()
			
			ret = ""
			for tag in tags:
				ret += tag + "\n"
			
			pluginHost.returnData(ret)

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
	ShowGeoTiffTags(args)
