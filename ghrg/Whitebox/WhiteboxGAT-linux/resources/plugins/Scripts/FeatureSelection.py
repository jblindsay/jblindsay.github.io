'''
 * Copyright (C) 2014 Dr. John Lindsay <jlindsay@uoguelph.ca>
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
from threading import Thread
from whitebox.ui.plugin_dialog import ScriptDialog
from java.awt.event import ActionListener
from whitebox.interfaces import MapLayer
from whitebox.interfaces.MapLayer import MapLayerType
from whitebox.geospatialfiles import ShapeFile
from whitebox.geospatialfiles import VectorLayerInfo
from whiteboxgis.user_interfaces import AttributesFileViewer

name = "FeatureSelection" 
descriptiveName = "Feature Selection" 
description = "This tool processes counts the number of records in a shapefile." 
toolboxes = ["VectorTools"] 
	
class FeatureSelection(ActionListener):
	def __init__(self, args, descriptiveName):
		if len(args) != 0:
			''' This will be followed if the parameters have already 
			been set, e.g. if the tool was called from another script. '''
			self.execute(args)
		else:
			''' The script dialog is created here. The first parameter 
			points the dialog to the Whitebox user interface. The second 
			parameter is the title of the dialog. In most cases this should 
			be set to the descriptive name of the tool. The third parameter
			is a listener for the OK button. When the user presses OK on 
			the dialog, actionPerformed method of this class will be informed.''' 
			self.sd = ScriptDialog(pluginHost, descriptiveName, self)	
			
			''' Sets the name of the helpfile, which should be an html file 
			contained in the Help directory. Ideally, it has the same name as 
			the tool's name.'''
			helpFile = self.__class__.__name__
			self.sd.setHelpFile(helpFile)
	
			''' This next line will enable the 'View Code' button on the dialog.'''
			self.sd.setSourceFile(os.path.abspath(__file__))
	
			''' Now add some components to the dialog. The following is a 
			description of some of the most commonly used dialog components.'''
	
			''' Use the addDialogFile method to add a File input/output component 
			to the dialog. The method has six parameters: 1) String description, 
			2) String labelText, 3) String dialogMode ("open" or "save"), 
			4) String filter (commonly "Raster Files (*.dep), DEP", 
			"Vector Files (*.shp), SHP", or "Whitebox Files (*.dep; *.shp), DEP, SHP"), 
			5) boolean showButton, 6) boolean makeOptional.'''
			self.sd.addDialogFile("Input vector layer", "Input Vector Layer:", "open", "Vector Files (*.shp), SHP", False, False)
			
			''' Resize the dialog to the standard size and display it '''
			self.sd.setSize(800, 400)
			self.sd.visible = True


	''' This method is the main part of your script. It is were the 
	main processing occurs and will be called on its own thread when 
	the user presses the OK button on the dialog (see actionPerformed
	below).'''
	def execute(self, args):
		try:
			''' Make sure that the args array has the expected number 
			of parameters. If the script is being called by another 
			script, this may not be the case, at which point, the script 
			should inform the user of the error and then end gracefully.'''
			if len(args) != 1:
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return

			inputfile = args[0]

			''' Feature selection is actually carried out in the attribute table dialog;
			this tool really just opens that dialog, feeding it the appropriate parameters. '''
			a = pluginHost.getAllMapLayers()

			found = False
			for k in a:
				if k.getLayerType() == MapLayerType.VECTOR:
						if k.getFileName() == inputfile:
							afv = AttributesFileViewer(pluginHost, False, k)
							afv.setActiveTab(3)
							afv.setSize(int(500 * 1.61803399), 500)
							afv.setVisible(True)
							found = True

			if not found:
				pluginHost.returnData(inputfile)

				a = pluginHost.getAllMapLayers()
				
				found = False
				for k in a:
					if k.getLayerType() == MapLayerType.VECTOR:
						if k.getFileName() == inputfile:
							afv = AttributesFileViewer(pluginHost, False, k)
							afv.setActiveTab(3)
							afv.setSize(int(500 * 1.61803399), 500)
							afv.setVisible(True)
							found = True
			
			if not found:
				pluginHost.showFeedback("There was an error while trying to open the Feature Selection dialog.")
				
		except Exception, e:
			print e
			pluginHost.showFeedback("Error during script execution.")
			''' alternatively, you many want to send the exception to 
			the pluginHost.logException() method '''
			return
		finally:
			pluginHost.updateProgress(0)


	''' This method is used as a listener for the dialog. If 
	the user presses the OK button, it will collect the 
	user-specified parameters then send them to the execute 
	method above. '''
	def actionPerformed(self, event):
		if event.getActionCommand() == "ok":
			args = self.sd.collectParameters()
			t = Thread(target=lambda: self.execute(args))
			t.start()

if args is None:
	''' Your script will have a variable called 'pluginHost' injected 
	into it. This variable is the Whitebox user interface and you 
	can use it to interact with the interface to do things like 
	show the user feedback, change the progress bar, check to see 
	if the user pressed the cancel button, display a file, etc.'''
	pluginHost.showFeedback("The arguments array has not been set.")
else:		
	FeatureSelection(args, descriptiveName)
