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

name = "LidarElevationSlice" 
descriptiveName = "LiDAR Elevation Slice" 
description = "Creates a new LAS file containing the points within an elevation range in an input LAS file" 
toolboxes = ["LidarTools"] 
	
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
			self.sd.addDialogFile("Input LAS file", "Input LAS File:", "open", "LAS Files (*.las), LAS", True, False)
			self.sd.addDialogFile("Output LAS file", "Output LAS File:", "save", "LAS Files (*.las), LAS", True, False)
			self.sd.addDialogDataInput("Minimum elevation in z-units", "Minimum elevation (Optional):", "", True, True)
			self.sd.addDialogDataInput("Maximum elevation in z-units", "Maximum elevation (Optional):", "", True, True)
			self.sd.addDialogDataInput("In-slice class value (0-31)", "In-slice Class Value (Optional):", "", True, True)
			self.sd.addDialogDataInput("Out-of-slice class value (0-31)", "Out-of-slice Class Value (Optional):", "", True, True)
			
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
			if len(args) != 6:
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return

			# read the input parameters
			inputfile = args[0]
			outputfile = args[1]
			min_z = float('-inf')
			if args[2].lower() != "not specified":
				min_z = float(args[2])
			max_z = float('inf')
			if args[3].lower() != "not specified":
				max_z = float(args[3])
			in_class = -1;
			if args[4].lower() != "not specified":
				in_class = int(args[4])
			out_class = -1;
			if args[5].lower() != "not specified":
				out_class = int(args[5])

			if min_z == float('-inf') and max_z == float('inf'):
				pluginHost.showFeedback("You must specify either a min. z, a max. z, or both.")
				return
			
			exe_path = pluginHost.getResourcesDirectory() + "plugins" + os.path.sep
			os.chdir(exe_path)

			(release, vendor, vminfo, osinfo) = jav()
			if "win" in osinfo[0].lower():
				ext = '.exe'
			else:
				ext = ''

			tool_name = "lidar_elevation_slice"
			cmd = "." + os.path.sep + "NativePlugins" + os.path.sep + "{}{}".format(tool_name, ext)
			cmd += ' -i=\"{}\"'.format(inputfile)
			cmd += ' -o=\"{}\"'.format(outputfile)
			cmd += ' -minz=\"{}\"'.format(min_z)
			cmd += ' -maxz=\"{}\"'.format(max_z)
			if in_class >= 0 and out_class >= 0:
				cmd += ' -class'
				cmd += ' -inclassval=\"{}\"'.format(in_class)
				cmd += ' -outclassval=\"{}\"'.format(out_class)
			cmd += ' -v'

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

			# display the output image
			pluginHost.returnData(outputfile)
			
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
