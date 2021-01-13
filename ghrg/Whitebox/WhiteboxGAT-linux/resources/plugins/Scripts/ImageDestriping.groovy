/*
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
def name = "ImageDestriping"
def descriptiveName = "Image De-striping"
def description = "Removes scan-line stripes in an image."
def toolboxes = ["ImageEnhancement"]

public class ImageDestriping implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public ImageDestriping(WhiteboxPluginHost pluginHost, 
		String[] args, def descriptiveName) {
		this.pluginHost = pluginHost
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
			def helpFile = "ImageDestriping"
			sd.setHelpFile(helpFile)
		
			// Specifying the source file allows the 'view code' 
			// button on the tool dialog to be displayed.
			def pathSep = File.separator
			def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + "ImageDestriping.groovy"
			sd.setSourceFile(scriptFile)
			
			// add some components to the dialog
			sd.addDialogFile("Input image file", "Input Image:", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Output file", "Output Raster File:", "save", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogDataInput("Threshold value.", "Threshold Value:", "5.0", true, false)
            sd.addDialogComboBox("Which direction are the stripes in the image?", "Striping Direction:", ["Horizontal", "Vertical"], 0)
			
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
			if (args.length != 4) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String inputFile = args[0]
			String outputFile = args[1]
			double diffThreshold = Double.parseDouble(args[2])
			String direction = "horizontal";
			if (args[3].toLowerCase().contains("vert")) {
				direction = "vertical"
			}

			// read the input image and PP vector files
			WhiteboxRaster image = new WhiteboxRaster(inputFile, "r")
			double nodata = image.getNoDataValue()
			int rows = image.getNumberRows()
			int cols = image.getNumberColumns()
			
			WhiteboxRaster output = new WhiteboxRaster(outputFile, "rw", 
  		  	  inputFile, DataType.FLOAT, nodata)

			double z, zN1, zN2, zOut, diff1, diff2
			int filterSize = 2
			int progress = 0
			int oldProgress = -1
			if (direction.equals("horizontal")) {
				for (int row = 0; row < rows; row++) {
	  				for (int col = 0; col < cols; col++) {
	  					zOut = image.getValue(row, col)
	  					if (zOut != nodata) {
	  						boolean flag = true
	  						for (int a = col - filterSize; a < col + filterSize; a++) {
	  							z = image.getValue(row, a)
	  							zN1 = image.getValue(row - 1, a)
	  							zN2 = image.getValue(row + 1, a)
	  							if (z != nodata && zN1 != nodata && zN2 != nodata) {
	  								diff1 = z - zN1
	  								diff2 = z - zN2
	  								if (Math.abs(diff1) < diffThreshold || 
	  								  Math.abs(diff2) < diffThreshold || 
	  								  (diff1 < 0 && diff2 >= 0) || 
	  								  (diff1 >= 0 && diff2 < 0)) {
	  									flag = false
	  									break
	  								}
	  							} else {
	  								flag = false
	  								break
	  							}
	  						}
	  						if (flag) {
	  							zN1 = image.getValue(row - 1, col)
	  							zN2 = image.getValue(row + 1, col)
	  							zOut = (zN1 + zN2) / 2.0
	  						}
	  					}
	  					
	  					output.setValue(row, col, zOut)
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
			} else { // direction = vertical
				for (int row = 0; row < rows; row++) {
	  				for (int col = 0; col < cols; col++) {
	  					zOut = image.getValue(row, col)
	  					if (zOut != nodata) {
	  						boolean flag = true
	  						for (int a = row - filterSize; a < row + filterSize; a++) {
	  							z = image.getValue(a, a)
	  							zN1 = image.getValue(a, col - 1)
	  							zN2 = image.getValue(a, col + 1)
	  							if (z != nodata && zN1 != nodata && zN2 != nodata) {
	  								diff1 = z - zN1
	  								diff2 = z - zN2
	  								if (Math.abs(diff1) < diffThreshold || 
	  								  Math.abs(diff2) < diffThreshold || 
	  								  (diff1 < 0 && diff2 >= 0) || 
	  								  (diff1 >= 0 && diff2 < 0)) {
	  									flag = false
	  									break
	  								}
	  							} else {
	  								flag = false
	  								break
	  							}
	  						}
	  						if (flag) {
	  							zN1 = image.getValue(row, col - 1)
	  							zN2 = image.getValue(row, col + 1)
	  							zOut = (zN1 + zN2) / 2.0
	  						}
	  					}
	  					
	  					output.setValue(row, col, zOut)
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
			}

			output.setDisplayMinimum(image.getDisplayMinimum())
			output.setDisplayMaximum(image.getDisplayMaximum())
			image.close()
			
			output.addMetadataEntry("Created by the " + descriptiveName + " tool.")
	        output.addMetadataEntry("Created on " + new Date())
			output.close()
	
			// display the output image
			pluginHost.returnData(outputFile)
	
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
	def f = new ImageDestriping(pluginHost, args, descriptiveName)
}
