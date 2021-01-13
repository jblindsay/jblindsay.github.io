# Copyright (C) 2017 Dr. John Lindsay <jlindsay@uoguelph.ca>
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import os
import sys
from sys import platform
from platform import java_ver as jav
import subprocess
import time
import math
from threading import Thread
from whitebox.ui.plugin_dialog import ScriptDialog
from java.awt.event import ActionListener
from whitebox.geospatialfiles import WhiteboxRaster
from whitebox.geospatialfiles.WhiteboxRasterBase import DataType

name = "ExportGeoTIFF" 
descriptiveName = "Export GeoTIFF" 
description = "Exports a Whitebox GAT Raster to GeoTIFF format (.tif)" 
toolboxes = ["IOTools"] 
	
class PluginTool(ActionListener):
	def __init__(self, args):
		if len(args) != 0:
			self.execute(args)
		else:
			''' Create a dialog for this tool to collect user-specified
			   tool parameters.''' 
			self.sd = ScriptDialog(pluginHost, descriptiveName, self)	
			
			''' Specifying the help file will display the html help
			// file in the help pane. This file should be be located 
			// in the help directory and have the same name as the 
			// class, with an html extension.'''
			self.sd.setHelpFile(name)
	
			''' Specifying the source file allows the 'view code' 
			// button on the tool dialog to be displayed.'''
			self.sd.setSourceFile(os.path.abspath(__file__))
	
			# add some components to the dialog '''
			self.sd.addDialogMultiFile("Select the input Whitebox raster files", "Input Whitebox Raster Files:", "Whitebox Raster (*.dep), DEP")
			
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
			if len(args) != 1:
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return

			# read the input parameters
			input_files = args[0][:-1]
			
			exe_path = pluginHost.getResourcesDirectory() + "plugins" + os.path.sep
			os.chdir(exe_path)

			(release, vendor, vminfo, osinfo) = jav()
			if "win" in osinfo[0].lower():
				ext = '.exe'
			else:
				ext = ''

			exe_name = "go-spatial"
			tool_name = "whitebox2geotiff"

			input_files_list = input_files.split(";")
			num_files = len(input_files_list)
			
			for input_file in input_files_list:
				cmd = "." + os.path.sep + "NativePlugins" + os.path.sep + "{0}{1}".format(exe_name, ext)
				cmd += ' -run=\"{}\"'.format(tool_name)
				cmd += ' -args=\"{0};{1}\"'.format(input_file, input_file.replace(".dep", ".tif"))

				ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, universal_newlines=True)
				
				while True:
					line = ps.stdout.readline()
					if line != '':
						if "%" in line:
							str_array = line.split(" ")
							label = line.replace(str_array[len(str_array)-1], "")
							progress = int(str_array[len(str_array)-1].replace("%", "").strip())
							pluginHost.updateProgress(label, progress)
						else:
							if "error" in line.lower():
								pluginHost.showFeedback("Error: {}".format(line))
							else:
								if not line.startswith("*"):
									pluginHost.updateProgress(line, 0)
					else:
						break

			pluginHost.showFeedback("Complete!")
		
		except Exception, e:
			print e
			pluginHost.showFeedback("An error has occurred during operation. See log file for details.")
			pluginHost.logException("Error in " + descriptiveName, e)
			return
		finally:
			# reset the progress bar
			pluginHost.updateProgress("Progress", 0)
			
if args is None:
	pluginHost.showFeedback("The arguments array has not been set.")
else:
	PluginTool(args)
