/*
 * Copyright (C) 2016 Dr. John Lindsay <jlindsay@uoguelph.ca>
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
def name = "MultiscaleElevationResidualIndex"
def descriptiveName = "Multiscale Elevation Residual Index"
def description = "Calculates the multi-scale elevation residual index for a DEM."
def toolboxes = ["ElevResiduals"]

public class MultiscaleElevationResidualIndex implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public MultiscaleElevationResidualIndex(WhiteboxPluginHost pluginHost, 
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
			def helpFile = "MultiscaleElevationResidualIndex"
			sd.setHelpFile(helpFile)
		
			// Specifying the source file allows the 'view code' 
			// button on the tool dialog to be displayed.
			def pathSep = File.separator
			def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + "MultiscaleElevationResidualIndex.groovy"
			sd.setSourceFile(scriptFile)
			
			// add some components to the dialog
			sd.addDialogFile("Input DEM file", "Input DEM:", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Output file", "Output Raster File:", "save", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogDataInput("Base value", "Base Value:", "1.5", true, false)
            sd.addDialogDataInput("Number of significant decimal places", "Significant Decimal Places:", "2", true, true)
            
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
			if (args.length < 3) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String inputFile = args[0]
			String outputFile = args[1]
			double baseValue = Double.parseDouble(args[2])
			if (baseValue <= 1) { baseValue = 1.01 }
			if (baseValue >= 2) { baseValue = 2.0 }
			int numSigDecimalPlaces = 2;
			if (args.length >= 4) {
				numSigDecimalPlaces = Integer.parseInt(args[3]);
				if (numSigDecimalPlaces < 0) { numSigDecimalPlaces = 0; }
				if (numSigDecimalPlaces > 8) { numSigDecimalPlaces = 8; }
			}
			double multiplier = Math.pow(10, numSigDecimalPlaces);
			
			int numScales
			int numCells
			int neighbourhoodSize
			
			// read the input image
			WhiteboxRaster image = new WhiteboxRaster(inputFile, "r")
			double nodata = image.getNoDataValue()
			int rows = image.getNumberRows()
			int cols = image.getNumberColumns()
			double minValue = image.getMinimumValue()
			double range = image.getMaximumValue() - minValue

			int rMax = (int)(Math.floor(Math.log(Math.min((double)(cols / 2.0), (double)(rows / 2.0))) / Math.log(baseValue)))
			int oldNeighbourhoodSize = -1
			ArrayList<Integer> radii = new ArrayList<>()
			for (int i = 0; i <= rMax; i++) {
				neighbourhoodSize = (int)(Math.floor(Math.pow(baseValue, i)))
				if (neighbourhoodSize > oldNeighbourhoodSize) {
					radii.add(neighbourhoodSize)
					oldNeighbourhoodSize = neighbourhoodSize
				}
			}

			numScales = radii.size()

			long[][] integralImage = new long[rows][cols]
			int[][] integralImageN = new int[rows][cols]
			
			// calculate the integral image
			int progress = 0
			int oldProgress = -1
			double z
			long sum, val
			int sumN
			for (int row = 0; row < rows; row++) {
				sum = 0
				sumN = 0
  				for (int col = 0; col < cols; col++) {
  					z = image.getValue(row, col)
  					if (z == nodata) {
  						z = 0
  					} else {
  						z = (z - minValue) * multiplier //(z - minValue) / range
  						sumN++
  					}
  					val = (long)(Math.round(z))
  					sum += val //z
  					if (row > 0) {
						integralImage[row][col] = sum + integralImage[row - 1][col]
						integralImageN[row][col] = (int)(sumN + integralImageN[row - 1][col])
  					} else {
						integralImage[row][col] = sum
						integralImageN[row][col] = (int)sumN
  					}
  				}
  				progress = (int)(100f * row / (rows - 1))
				if (progress > oldProgress) {
					pluginHost.updateProgress("Loop 1 of 2:", progress)
					oldProgress = progress
				}
				// check to see if the user has requested a cancellation
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled")
					return
				}
			}

			WhiteboxRaster output = new WhiteboxRaster(outputFile, "rw", 
  		  	  inputFile, DataType.FLOAT, 0.0)
			output.setPreferredPalette("grey.plt")
			
			int x1, x2, y1, y2
			double a, b, c, d, outValue
			double value

			oldProgress = -1
			for (int row = 0; row < rows; row++) {
				for (int col = 0; col < cols; col++) {
					z = image.getValue(row, col)
  					if (z != nodata) {
  						neighbourhoodSize = 1
  						value = 0
						for (int r : radii) {
							y1 = row - r
							if (y1 < 0) { y1 = 0 }
							if (y1 >= rows) { y1 = rows - 1 }
			
							y2 = row + r
							if (y2 < 0) { y2 = 0 }
							if (y2 >= rows) { y2 = rows - 1 }
						
		  					x1 = col - r
							if (x1 < 0) { x1 = 0 }
							if (x1 >= cols) { x1 = cols - 1 }
		
							x2 = col + r
							if (x2 < 0) { x2 = 0 }
							if (x2 >= cols) { x2 = cols - 1 }

							a = integralImage[y1][x1] / multiplier
							b = integralImage[y1][x2] / multiplier
							c = integralImage[y2][x2] / multiplier
							d = integralImage[y2][x1] / multiplier

							numCells = (int)(integralImageN[y2][x2] + integralImageN[y1][x1] - integralImageN[y1][x2] - integralImageN[y2][x1])
	
							//outValue = z - ((c + a - b - d) / numCells * range + minValue)
							outValue = z - ((c + a - b - d) / numCells + minValue)
							if (outValue > 0) {
								value++
							}
						}
						output.setValue(row, col, value / numScales)
  					} else {
  						output.setValue(row, col, nodata)
  					}
  				}
  				progress = (int)(100f * row / (rows - 1))
				if (progress > oldProgress) {
					pluginHost.updateProgress("Loop 2 of 2:", progress)
					oldProgress = progress
				}
				// check to see if the user has requested a cancellation
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled")
					return
				}
			}

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
	def f = new MultiscaleElevationResidualIndex(pluginHost, args, descriptiveName)
}
