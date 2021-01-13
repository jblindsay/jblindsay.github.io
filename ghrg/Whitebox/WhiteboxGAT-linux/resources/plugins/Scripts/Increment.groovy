/*
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
 */

import java.awt.event.ActionListener
import java.awt.event.ActionEvent
import java.util.concurrent.Future
import java.util.concurrent.*
import java.util.Date
import java.util.ArrayList
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterInfo
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "Increment"
def descriptiveName = "Increment"
def description = "Increments the values in a raster grid by the values in a second grid or a constant."
def toolboxes = ["MathTools"]

public class Increment implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd
	private String descriptiveName
	private String name
	
	public Increment(WhiteboxPluginHost pluginHost, 
		String[] args, def name, def descriptiveName) {
		this.pluginHost = pluginHost
		this.name = name
		this.descriptiveName = descriptiveName
			
		if (args.length > 0) {
			execute(args)
		} else {
			// Create a dialog for this tool to collect user-specified
			// tool parameters.
			sd = new ScriptDialog(pluginHost, descriptiveName, this)	
		
			// Specifying the help file will display the html help
			// file in the help pane. This file should be be located 
			// in the help directory and have the same name as the 
			// class, with an html extension.
			sd.setHelpFile(name)
		
			// Specifying the source file allows the 'view code' 
			// button on the tool dialog to be displayed.
			def pathSep = File.separator
			def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + name + ".groovy"
			sd.setSourceFile(scriptFile)
			
			// add some components to the dialog
			sd.addDialogFile("Input image file", "Image To Be Incremented:", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Increment value source (raster or constant value)", "Increment Value (Raster File Or Constant Value):", "open", "Raster Files (*.dep), DEP", true, false)
			
			// resize the dialog to the standard size and display it
			sd.setSize(800, 400)
			sd.visible = true
		}
	}

	// The CompileStatic annotation can be used to significantly
	// improve the performance of a Groovy script to nearly 
	// that of native Java code.
	@CompileStatic
	private void execute(String[] args) {
		try {
			double z, incrementVal
			int progress, oldProgress
			if (args.length != 2) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String inputFile = args[0]
			def file = new File(args[1])
			boolean incrementValsConstant = !(file).exists()
			String incrementFile = ""
			if (incrementValsConstant) {
				incrementVal = Double.parseDouble(file.getName().replace(".dep", ""))
			} else {
				incrementFile = args[1]
			}
			
			// read the input image and PP vector files
			WhiteboxRaster image = new WhiteboxRaster(inputFile, "rw")
			double nodata = image.getNoDataValue()
			int rows = image.getNumberRows()
			int cols = image.getNumberColumns()

			oldProgress = -1
			if (incrementValsConstant) {
				for (int row = 0; row < rows; row++) {
	  				for (int col = 0; col < cols; col++) {
	  					z = image.getValue(row, col)
	  					if (z != nodata) {
							image.incrementValue(row, col, incrementVal)
	  					}
	  				}
	  				progress = (int)(100f * row / rows)
					if (progress > oldProgress) {
						pluginHost.updateProgress(progress)
						oldProgress = progress
					}
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			} else {
				WhiteboxRaster incrementImage = new WhiteboxRaster(incrementFile, "r")
				double incrementNodata = incrementImage.getNoDataValue()
				if (incrementImage.getNumberRows() != rows || incrementImage.getNumberColumns() != cols)  {
					pluginHost.showFeedback("The input image and increment image do not have the same dimensions.")
					return
				}
				
				for (int row = 0; row < rows; row++) {
	  				for (int col = 0; col < cols; col++) {
	  					z = image.getValue(row, col)
	  					if (z != nodata) {
	  						incrementVal = incrementImage.getValue(row, col)
	  						if (incrementVal != incrementNodata) {
								image.incrementValue(row, col, incrementVal)
	  						}
	  					}
	  				}
	  				progress = (int)(100f * row / rows)
					if (progress > oldProgress) {
						pluginHost.updateProgress(progress)
						oldProgress = progress
					}
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
				incrementImage.close()
			}
			image.flush()
			image.findMinAndMaxVals()
			image.resetDisplayMinMaxValues()
			image.deleteStatsFile()
			image.close()
			
			// display the image
			pluginHost.returnData(inputFile)
	
		} catch (OutOfMemoryError oe) {
            pluginHost.showFeedback("An out-of-memory error has occurred during operation.")
	    } catch (Exception e) {
	        pluginHost.showFeedback("An error has occurred during operation. See log file for details.")
	        pluginHost.logException("Error in " + descriptiveName, e)
        } finally {
        	// reset the progress bar
        	pluginHost.updateProgress(0)
        }
	}
	
	@Override
    public void actionPerformed(ActionEvent event) {
    	if (event.getActionCommand().equals("ok")) {
    		final def args = sd.collectParameters()
			sd.dispose()
			final Runnable r = new Runnable() {
            	@Override
            	public void run() {
                	execute(args)
            	}
        	}
        	final Thread t = new Thread(r)
        	t.start()
    	}
    }
}

if (args == null) {
	pluginHost.showFeedback("Plugin arguments not set.")
} else {
	def f = new Increment(pluginHost, args, name, descriptiveName)
}
