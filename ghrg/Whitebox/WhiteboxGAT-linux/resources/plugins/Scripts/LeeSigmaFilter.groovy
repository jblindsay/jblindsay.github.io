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
def name = "LeeSigmaFilter"
def descriptiveName = "Lee (Sigma) Filter"
def description = "Performs a Lee (Sigma) smoothing filter on an image."
def toolboxes = ["Filters"]

public class LeeSigmaFilter implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public LeeSigmaFilter(WhiteboxPluginHost pluginHost, 
		String[] args, def name, def descriptiveName) {
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
			sd.setHelpFile(name)
		
			// Specifying the source file allows the 'view code' 
			// button on the tool dialog to be displayed.
			def pathSep = File.separator
			def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + name + ".groovy"
			sd.setSourceFile(scriptFile)
			
			// add some components to the dialog
			sd.addDialogFile("Input raster image file", "Input Raster:", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Output file", "Output Raster File:", "save", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogDataInput("Neighbourhood size in X dimension", "Neighbourhood Size in X (cells):", "5", true, false)
            sd.addDialogDataInput("Neighbourhood size in Y dimension", "Neighbourhood Size in Y (cells):", "5", true, false)
            sd.addDialogDataInput("Sigma", "Sigma:", "", true, false)
            sd.addDialogDataInput("M-threshold", "M-threshold:", "", true, false)
            sd.addDialogCheckBox("Reflect values at image edges?", "Reflect values at image edges?", true)
			
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
			if (args.length != 7) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String inputFile = args[0]
			String outputFile = args[1]
			int filterSizeX = Integer.parseInt(args[2])
			if (filterSizeX < 3) { filterSizeX = 3 }
			int filterSizeY = Integer.parseInt(args[3])
			if (filterSizeY < 3) { filterSizeY = 3 }
			double sigma = Double.parseDouble(args[4])
			sigma = sigma * 2
			int M = Integer.parseInt(args[5])
			if (M < 0) { M = 0 }
			boolean reflectAtEdge = Boolean.parseBoolean(args[6])

			//the filter dimensions must be odd numbers such that there is a middle pixel
            if (Math.floor(filterSizeX / 2d) == (filterSizeX / 2d)) {
                pluginHost.showFeedback("Filter dimensions must be odd numbers. The specified filter x-dimension" + 
                        " has been modified.");
                
                filterSizeX++;
            }
            if (Math.floor(filterSizeY / 2d) == (filterSizeY / 2d)) {
                pluginHost.showFeedback("Filter dimensions must be odd numbers. The specified filter y-dimension" + 
                        " has been modified.");
                filterSizeY++;
            }
			
			//fill the filter kernel cell offset values
            int numPixelsInFilter = filterSizeX * filterSizeY;
            int[] dX = new int[numPixelsInFilter];
            int[] dY = new int[numPixelsInFilter];

            int midPointX = (int)Math.floor(filterSizeX / 2);
            int midPointY = (int)Math.floor(filterSizeY / 2);

            int a = 0;
            for (int row = 0; row < filterSizeY; row++) {
                for (int col = 0; col < filterSizeX; col++) {
                    dX[a] = col - midPointX;
                    dY[a] = row - midPointY;
                    a++;
                 }
            }

			// these are the filter kernel cell offsets for the
			// 3 x 3 neighbourhood that is used to calculate
			// the average of the immediate neighbourhood when n < M
            int[] dX1 = [ 1, 1, 1, 0, -1, -1, -1, 0 ]
			int[] dY1 = [ -1, 0, 1, 1, 1, 0, -1, -1 ]
			
			
			// read the input image file
			WhiteboxRaster image = new WhiteboxRaster(inputFile, "r")
			double nodata = image.getNoDataValue()
			int rows = image.getNumberRows()
			int cols = image.getNumberColumns()
			image.isReflectedAtEdges = reflectAtEdge

			// create the output image
			WhiteboxRaster output = new WhiteboxRaster(outputFile, "rw", 
  		  	  inputFile, DataType.FLOAT, nodata)
  		  	//output.setNoDataValue(nodata)
			output.setPreferredPalette(image.getPreferredPalette())
			
			
			int progress, oldProgress, n
			double z, zN, sum, upperValue, lowerValue
			for (int row = 0; row < rows; row++) {
				for (int col = 0; col < cols; col++) {
  					z = image.getValue(row, col)
  					upperValue = z + sigma
  					lowerValue = z - sigma
  					if (z != nodata) {
  						n = 0
                        sum = 0
                        for (a = 0; a < numPixelsInFilter; a++) {
                            zN = image.getValue(row + dY[a], col + dX[a])
                           	if (zN >= lowerValue && zN <= upperValue && zN != nodata) {
                                n++
                                sum += zN
                            }
                        }
                        
                        if (n > M) {
                            output.setValue(row, col, sum / n)
                        } else {
                        	n = 0
                        	sum = 0
                        	for (a = 0; a < 8; a++) {
                        		zN = image.getValue(row + dY1[a], col + dX1[a])
	                            if (zN != nodata) {
	                                n++
	                                sum += zN
	                            }
                        	}
                        	if (n > 0) {
                            	output.setValue(row, col, sum / n)
                        	}
                        }       
  					}
  				}
  				progress = (int)(100f * row / (rows - 1))
				if (progress != oldProgress) {
					pluginHost.updateProgress(progress)
					oldProgress = progress
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}

			
			image.close()

			output.flush()
			output.findMinAndMaxVals()
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
	def f = new LeeSigmaFilter(pluginHost, args, name, descriptiveName)
}
